#include <plan_manage/nmpc_utils.h>

constexpr double PI = 3.1415926;

namespace resilient_planner
{

    /// 初始化 避障约束的参数
  NMPCSolver::NMPCSolver()
  {
    //construct model matrices
    At_.setZero();
    Bt_.setZero();
    Dt_.setZero();

    // some pre initializations of the matrix
    At_(0, 3) = 1.0;
    At_(1, 4) = 1.0;
    At_(2, 5) = 1.0;

    Bt_(6, 0) = 1.0;
    Bt_(7, 1) = 1.0;
    Bt_(8, 2) = 1.0;

    Dt_(3, 0) = 1.0;
    Dt_(4, 1) = 1.0;
    Dt_(5, 2) = 1.0;

    Kt_ <<  -2.0,  5.0,  0.0, -1.0,  4.0,  0.0, -8.0,  0.0,  0.0,
            -5.0, -2.0,  0.0, -4.0, -1.0,  0.0,  0.0, -8.0,  0.0,
            -2.0, -2.0,  0.0, -1.0, -1.0,  0.0,  0.0,  0.0, -8.0,
             0.0,  0.0, -8.0,  0.0,  0.0, -6.0,  0.0,  0.0,  0.0;

    external_acc_.setZero();
    end_v_.setZero();

    cmd_status_ = CMD_STATUS::INIT_POSITION;
  }

  /// 初始化ROS空间: 包括前端的路径搜索、后端优化参数
  void NMPCSolver::initROS(ros::NodeHandle &nh, OccMap::Ptr &env_ptr_)
  {

    /*  kino a* intial  */
    kino_path_finder_.reset(new KinodynamicAstar);
    kino_path_finder_->setParam(nh);
    kino_path_finder_->init();
    kino_path_finder_->intialGridMap(env_ptr_);

    ellipsoid_pub_ = nh.advertise<decomp_ros_msgs::EllipsoidArray>("ellipsoid_array", 1, true);  ///
    poly_pub_ = nh.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedron_array", 1, true);  ///

    ref_marker_pub_ = nh.advertise<visualization_msgs::Marker>("ref_traj", 1, true);  /// ref traj
    nmpc_marker_pub_ = nh.advertise<visualization_msgs::Marker>("nmpc_traj", 1, true); /// nmpc traj
    cmd_vis_pub_ = nh.advertise<visualization_msgs::Marker>("position_cmd_vis", 10, true); /// pos_cmd

    cloud_sub_ = nh.subscribe("/occ_map/local_view_cloud", 1, &NMPCSolver::cloudCallback, this);   /// 局部点云数据

    double ego_r, ego_h, drag_coeff, ext_noise_bound;

    /* weight */
    double w_stage_wp, w_stage_input, w_terminal_wp, w_terminal_input,w_input_rate,
           w_final_terminal_wp, w_final_terminal_input,w_final_stage_wp, w_final_stage_input;

    nh.param("nmpc/w_stage_wp", w_stage_wp, 15.0);
    nh.param("nmpc/w_stage_input", w_stage_input, 3.0);
    nh.param("nmpc/w_terminal_wp", w_terminal_wp, 15.0);
    nh.param("nmpc/w_terminal_input", w_terminal_input, 0.0);
    nh.param("nmpc/w_input_rate", w_input_rate, 80.0);
    nh.param("nmpc/w_final_stage_wp", w_final_stage_wp, 20.0);
    nh.param("nmpc/w_final_stage_input", w_final_stage_input, 5.0);
    nh.param("nmpc/w_final_terminal_wp", w_final_terminal_wp, 20.0);
    nh.param("nmpc/w_final_terminal_input", w_final_terminal_input, 5.0);

    nh.param("nmpc/ego_r", ego_r, 0.27);
    nh.param("nmpc/ego_h", ego_h, 0.0425);
    nh.param("nmpc/g_acc", g, 9.81);
    nh.param("nmpc/mass", mass, 0.74);
    nh.param("nmpc/drag_coefficient", drag_coeff, 0.33);
    nh.param("nmpc/ext_noise_bound", ext_noise_bound, 0.5);  ///

    nh.param("search/max_tau", max_tau_, 0.5);

    ROS_INFO_STREAM("ego_r: " << ego_r);
    ROS_INFO_STREAM("ego_h: " << ego_h);
    ROS_INFO_STREAM("max_tau: " << max_tau_);
    ROS_INFO_STREAM("mass: " << mass);

    // use rotors simulator
    cmd_timer_ = nh.createTimer(ros::Duration(0.01), &NMPCSolver::cmdTrajCallback, this); // 5ms
    traj_cmd_pub_ = nh.advertise<trajectory_msgs::MultiDOFJointTrajectory>("/pos_cmd", 50);

    // ego info
    ego_size_ << ego_r * ego_r, 0.0, 0.0,
        0.0, ego_r * ego_r, 0.0,
        0.0, 0.0, ego_h * ego_h;

    //drag_coefficient_matrix
    /// 没有考虑z轴
    drag_coefficient_matrix_ << drag_coeff, 0.0, 0.0,
        0.0, drag_coeff, 0.0,
        0.0, 0.0, 0.0;

    /// disturbance vector
    /// 干扰向量
    w_ << ext_noise_bound, ext_noise_bound, ext_noise_bound;

    nmpc_forces_solver_normal_.setParasNormal(w_stage_wp, w_stage_input, w_input_rate, w_terminal_wp, w_terminal_input);  /// normal:正常
    nmpc_forces_solver_final_.setParasFinal(w_final_stage_wp, w_final_stage_input, w_input_rate,  /// 最后一段轨迹
                                            w_final_terminal_wp, w_final_terminal_input);

    ROS_INFO_STREAM("[NMPCSolver]Finish the intialization of the nmpc solver! \n");
  }


  /// 获得kino搜索出的某个路径点
  void NMPCSolver::getCurTraj(int index)
  {
    double index_time = index * Ts_ + (mpc_start_time_ - kino_start_time_).toSec();
    unsigned int kino_index = (int)(index_time / Ts_);

    /// decide the reference point
    if (kino_index + 1 < kino_size_)
    {
      ref_pos_ = kino_path_[kino_index] + fmod(index_time, Ts_) / Ts_ * (kino_path_[kino_index + 1] - kino_path_[kino_index]);
    }
    else
    {
      ref_pos_ = kino_path_[kino_size_ - 1];
    }

    /// decide the forward point
    if (kino_index + 5 < kino_size_)
    {
      forward_pos_ = kino_path_[kino_index + 5];
    }
    else
    {
      forward_pos_ = kino_path_[kino_size_ - 1];
    }

    calculate_yaw(ref_pos_, forward_pos_);  /// 计算偏航角

    if (index == 0 && (ref_pos_ - mpc_output_.at(1).segment(8, 3)).norm() > 1.0)
    {
      ROS_INFO_STREAM("Hard to follow the reference ! replan !\n");
      kino_replan_ = true;
    }
    
  }

  /// use kinodynamic a* to generate a path
  bool NMPCSolver::getKinoPath(Eigen::VectorXd &stateOdom, Eigen::Vector3d end_pt,
                               Eigen::Vector3d external_acc, bool replan)
  {
    end_pt_ = end_pt;
    external_acc_ = external_acc;

    start_pt_ = stateOdom.segment(0, 3);
    start_v_ = stateOdom.segment(3, 3);
    start_a_.setZero();

    kino_path_finder_->reset();
    kino_path_finder_->updateExternalAcc(external_acc_);

    int status;

    if (replan && exit_code == 1)  /// 需要重规划 且上一次规划成功
    { // need replan and mpc success

       /// 如果上一次规划成功，就用上一次的规划成功的信息
      double t_cur = (ros::Time::now() - pre_mpc_start_time_).toSec();
      int cur_index = (int)(t_cur / Ts_);
      Eigen::VectorXd mpcq;

      if (cur_index < planning_horizon_ - 1 && t_cur >= 0.0)
      {
        mpcq = pre_mpc_output_.at(cur_index) + fmod(t_cur, Ts_) / Ts_ * (pre_mpc_output_.at(cur_index + 1) - pre_mpc_output_.at(cur_index));

        Eigen::Vector3d temp_start_pt, temp_start_v, temp_start_a;

        temp_start_pt << mpcq(8), mpcq(9), mpcq(10);
        temp_start_v << mpcq(11), mpcq(12), mpcq(13);

        Eigen::Vector3d odom_euler(mpcq(14), mpcq(15), mpcq(16));
        Eigen::Vector3d thrust_b(0.0, 0.0, mpcq(3));
        Eigen::Vector3d thrustn_w = eulerToRot(odom_euler) * thrust_b / mass;

        temp_start_a << thrustn_w(0), thrustn_w(1), thrustn_w(2) - g;

        status = kino_path_finder_->search(temp_start_pt, temp_start_v, temp_start_a, end_pt_, end_v_, true);
      }
      else
      {
        status = kino_path_finder_->search(start_pt_, start_v_, start_a_, end_pt_, end_v_, true);
      }
    }
    else  /// 如果上一次没有规划成功，就用当前里程计的信息
    {
      status = kino_path_finder_->search(start_pt_, start_v_, start_a_, end_pt_, end_v_, true);
    }

    if (status == KinodynamicAstar::NO_PATH)
    {
      std::cout << "[kino replan]: kinodynamic search fail!" << std::endl;

      // retry searching with discontinuous initial state
      kino_path_finder_->reset();
      status = kino_path_finder_->search(start_pt_, start_v_, start_a_, end_pt, end_v_, false);
      if (status == KinodynamicAstar::NO_PATH)
      {
        std::cout << "[kino replan]: Can't find path." << std::endl;
        return false;
      }
      else
      {
        std::cout << "[kino replan]: retry search success." << std::endl;
      }
    }
    else
    {
      std::cout << "[kino replan]: kinodynamic search success." << std::endl;
    }

    kino_path_ = kino_path_finder_->getKinoTraj(Ts_);  ///
    kino_size_ = kino_path_.size();

    switch_to_final = false;
    kino_start_time_ = ros::Time::now();
    ROS_INFO_STREAM("kino_start_time_ : \n"
                    << kino_start_time_);
    cmd_status_ = CMD_STATUS::PUB_TRAJ;   /// kino搜索成功 就发布nmpc的信息
    pub_end_ = false;  ///

    return true;
  }



  /// 初始化偏航角  定义init_yaw_dot_的数值和cmd_status_变量的数值
  void NMPCSolver::callInitYaw(Eigen::VectorXd odom, double init_yaw)
  {
    ROS_INFO_STREAM("[NMPCSolver] callInitYaw");
    realOdom_ = odom;
    init_yaw_ = init_yaw;
    pub_end_ = false;
    change_yaw_time_ = ros::Time::now();

    init_yaw_dot_ = init_yaw_ - realOdom_(8);  ///

    if (init_yaw_dot_ > PI)
    {
      //init_yaw_dot_ = 2 * PI - init_yaw_dot_;     /// 有问题
        init_yaw_dot_ = - 2 * PI + init_yaw_dot_;
    }
    else if (init_yaw_dot_ < -PI)
    {
        init_yaw_dot_ = init_yaw_dot_ + 2 * PI;
    }

    double max_yaw_dot = 0.4 * PI;  /// 限制了偏航角速度 没有限制偏航角最大数

    if (init_yaw_dot_ > max_yaw_dot)
    {
      init_yaw_dot_ = max_yaw_dot;
    }
    else if (init_yaw_dot_ < -max_yaw_dot)
    {
      init_yaw_dot_ = -max_yaw_dot;
    }

    ROS_INFO_STREAM("init_yaw_dot_ : \n"
                    << init_yaw_dot_);
    cmd_status_ = CMD_STATUS::ROTATE_YAW;   ///
  }


  ///
  void NMPCSolver::initMPCOutput()
  {
    mpc_output_.clear();

    Eigen::VectorXd mpc_row(var);
    // control input is
    // state is: position, velocity and euler angle
    mpc_row << 0.0, 0.0, 0.0, real_thrust_c_,
        0.0, 0.0, 0.0, real_thrust_c_,
        stateMpc_(0), stateMpc_(1), stateMpc_(2), // position
        stateMpc_(3), stateMpc_(4), stateMpc_(5), // velocity
        stateMpc_(6), stateMpc_(7), stateMpc_(8); // euler angle

    for (int i = 0; i < planning_horizon_ + 1; i++)
    {
      mpc_output_.push_back(mpc_row);   /// 初始化的时候每一步都一样
    }

    initialized_output_ = true;  /// 初始化完成的标志
    pre_mpc_output_ = mpc_output_; ///
    ROS_INFO_STREAM("[NMPCSolver::initMPCOutput] Finish the intialization of mpc output! ");
  }

  /// 验证避障约束  E：Q_ext^(1/2)
  int NMPCSolver::getSikangConst(Eigen::Matrix3d E)
  {
    int index = poly_constraints_.size();
    if (index >= 1)
    {
      /// directly check the reference ellipsoid shape within the polytope
      /// 直接检查是否满足约束
      auto temp_A = (poly_constraints_.back()).A();
      auto temp_b = (poly_constraints_.back()).b();
      bool flag = true;

      for (unsigned int j = 0; j < temp_b.size(); j++)
      {
        Eigen::MatrixXd temp_b_addition = E * (temp_A.row(j)).transpose();
        // with little inflation
        if (temp_A.row(j) * ref_pos_ - (temp_b(j) - 1.1 * temp_b_addition.norm()) > 0)
        {
          flag = false;
          break;
        }
      }

      // the point in the polytope
      if (flag)
        return index - 1; /// same as last time poly  都满足约束 那么和上一次的多面体不用更新
    }

    /// 不是全部满足的情况
    vec_Vec3f seed_path;
    seed_path.push_back(ref_pos_);

    Vec3f seed_point2;
    ///  注意怎么处理的
    seed_point2 << ref_pos_[0] + 0.1*cos(ref_yaw_), ref_pos_[1] + 0.1*sin(ref_yaw_), ref_pos_[2];  /// 第二参考点
    seed_path.push_back(seed_point2);

    EllipsoidDecomp3D decomp_util;  /// 第三方库
    decomp_util.set_obs(vec_obs_);    /// 输入的是点云向量
    decomp_util.set_local_bbox(Vec3f(2, 2, 1));  /// 局部边界框
    decomp_util.dilate(seed_path);  ///
    auto seed_poly = decomp_util.get_polyhedrons();  /// 获取多面体

    vec_E<LinearConstraint3D> css = decomp_util.get_constraints();  /// 获得约束
    polys.push_back(seed_poly[0]);  /// 增加一个多面体
    poly_constraints_.push_back(css[0]);  /// 增加一个多面体约束

    return index;
  }

  ///
  void NMPCSolver::displayPoly()
  {
    decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(polys);
    poly_msg.header.frame_id = "world";
    poly_pub_.publish(poly_msg);
  }

  ///
  void NMPCSolver::getKinoTraj(std::vector<Eigen::Vector3d> &kino_path)
  {
    kino_path = kino_path_;
  }

  //  1 --- success
  //  0 --- at global end
  // -1 --- reaching global end process but not finished
  // -2 --- need replan
  // -3 --- odom far away from predict state
  ///  求解 NMPC问题
  int NMPCSolver::solveNMPC(Eigen::VectorXd &stateMpc, Eigen::Vector3d external_acc)
  {
    external_acc_ = external_acc; // update external forces

    if (cmd_status_ == CMD_STATUS::WAIT) {
        return 0;
    }
    else if (pub_end_) {
        return -1;
    }

    mpc_start_time_ = ros::Time::now();
    pre_mpc_start_time_ = mpc_start_time_;
    stateMpc_ = stateMpc;

    // update the mpc output
    if (!initialized_output_ || exit_code != 1)
      initMPCOutput();  ///

    // for visualization
    poly_constraints_.clear();
    polys.clear();
    ellipsoid_matrices_.clear();
    ref_total_pos_.clear();
    ref_total_yaw_.clear();

    /*** intialization !!!  ***/
    clock_t start, finish1, finish2;
    double totalTime = 0;
    start = clock();
    bool normal_finish = false, update_result = false;

    Eigen::VectorXd poly_indices(planning_horizon_);
    //==============================================================================================//
    /// 获得 ref_total_pos_ ref_total_yaw_  ellipsoid_matrices_  poly_constraints_  poly_indices的数值 ，用于优化
    setFORCESParams(poly_indices);   /// 很重要！！！
    //==============================================================================================//
    if (!switch_to_final)
    {
      std::cout << "[NMPC NORMAL] Set up the solver ... " << std::endl;
      exit_code = nmpc_forces_solver_normal_.solveNormal(mpc_output_, external_acc_,
                                                         ref_total_pos_, ref_total_yaw_,
                                                         ellipsoid_matrices_, poly_constraints_, poly_indices);
      normal_finish = true;
    }
    else
    {
      std::cout << "[NMPC FINAL] Set up the solver ... " << std::endl;
      exit_code = nmpc_forces_solver_final_.solveFinal(mpc_output_, external_acc_,
                                                       ref_total_pos_, ref_total_yaw_,
                                                       ellipsoid_matrices_, poly_constraints_, poly_indices);
    }
    //============================================================================================================//
    /* check whether or not to use the results*/
    if (exit_code == 1) /// success  求解成功
    {
      // reset the variables//
      fail_count_ = 0;
      replan_count_ = 0;
      update_result = true;  /// 更新结果
    }
    else
    {
      fail_count_ += 1;
      if (replan_count_ > 3 && exit_code == 0)  /// 求解超时且kino重规划次数大于3（kino-opt-opt-opt-kin0-opt-opt-opt-kino-opt-opt-opt-kino-opt-opt-opt-kino-opt_replan_count=4 满足要求）
      {
        fail_count_ = 0;
        replan_count_ = 0;
        update_result = true;
      }
      else if (fail_count_ > 2)  /// 失败的次数大于2  kino重搜索
      {
        fail_count_ = 0;
        replan_count_ += 1;
        ROS_WARN("[NMPCsolver] MPC Fails too much time, replan !");
        kino_replan_ = true;  /// kino重规划
      }
    }

    if (update_result)
    {
      if (normal_finish){
          nmpc_forces_solver_normal_.updateNormal(mpc_output_);  /// 更新结果于mpc_output_
      }
      else {
          nmpc_forces_solver_final_.updateFinal(mpc_output_);  ///
      }

      updateFORCESResults(); ///
    }

    finish1 = clock();  /// clock()是C/C++中的计时函数，而与其相关的数据类型是clock_t
    /// CLOCKS_PER_SEC是标准c的time.h头函数中宏定义的一个常数，表示一秒钟内CPU运行的时钟周期数，用于将clock()函数的结果转化为以秒为单位的量，但是这个量的具体值是与操作系统相关的
    totalTime = (double)(finish1 - start) / CLOCKS_PER_SEC * 1000;
    printf("\033[34m[NMPCsolver] total time(ms)=%5.3f\n\033[0m", totalTime);
    //=========================================================================================================================//
    /// 值得借鉴
    /*    check replan or reach the global end    */
    Eigen::Vector3d ref_end = (mpc_output_.at(19)).segment(8, 3);  /// 第 planning_horizon-1组数据中的位置
    int max_index = (int)((planning_horizon_ * Ts_ + (mpc_start_time_ - kino_start_time_).toSec()) / Ts_);

    /// 搜索出的路径距离目标还比较远，同时nmpc优化出的路径点数目大于搜索出路径点数一半，很可能
    /// 下一次搜索失败（到达局部目标点，或者说局部目标点没有更新）
    if (max_index > 0.5 * kino_size_ && (end_pt_ - kino_path_[kino_size_ - 1]).norm() > 0.7)
    {
      ROS_WARN("[NMPCsolver] Reach the local end, replan!");
      kino_replan_ = true;
    }

    /* check whether to switch mode  */
    /// 接近最终的目标点
    if (max_index >= kino_size_ || (ref_end - end_pt_).norm() < 1.0)
      switch_to_final = true;
    // --------------------------check return module ----------------------//

    /// 优化的数据和里程计的数据差的太远   或许是外力太大无人机实际位置在受力瞬间产生较大的偏移
    /* check for odom and predict state  */
    if (((mpc_output_.at(1)).segment(8, 3) - stateMpc_.segment(0, 3)).norm() > 2.0)
    {
      std::cout << "[NMPCsolver] stateMpc_.segment(0, 3) ：" << stateMpc_.segment(0, 3) << std::endl;

      std::cout << "[NMPCsolver] (mpc_output_.at(1)).segment(8, 3)：" << (mpc_output_.at(1)).segment(8, 3)<< std::endl;
           
      cmd_status_ = CMD_STATUS::WAIT;
      return -3;
    }

    /// 与设定的目标点还要一段距离
    /* check whether reach the global end  */
    if ((ref_end - end_pt_).norm() < 0.15)
    {
      std::cout << "[NMPCsolver] Publish end position ：" << end_pt_ << std::endl;
      pub_end_ = true;
      cmd_end_pt_ = end_pt_;
      return -1;
    }

    /* check whether to replan  */
    if (kino_replan_)
    {
      kino_replan_ = false;
      return -2;
    }

    return 1; //SUCCESS
  }

  ///  两个矩阵：四旋翼本身Q1+外界扰动Q2
  void NMPCSolver::setFORCESParams(Eigen::VectorXd &poly_indices)
  {
    last_yaw_ = mpc_output_.at(1)(16);
    Eigen::MatrixXd Q_init = pow(epsilon, 2) * Eigen::MatrixXd::Identity(nx, nx);  /// 初始化安全椭球边界
    Eigen::Matrix3d Q, Q1, Q2;

    for (int i = 0; i < planning_horizon_; i++)
    {
      getCurTraj(i);  //
      ref_total_pos_.push_back(ref_pos_);
      ref_total_yaw_.push_back(ref_yaw_);

      /// use predicted euler angle
      Eigen::Vector3d euler_cur(mpc_output_.at(i)(14), mpc_output_.at(i)(15), mpc_output_.at(i)(16));
      Eigen::Vector3d vel_cur(mpc_output_.at(i)(11), mpc_output_.at(i)(12), mpc_output_.at(i)(13));

      Eigen::Matrix3d R_cur = updateMatrix(euler_cur, vel_cur, mpc_output_.at(i)(3));  /// 更新旋转矩阵 系统矩阵以及其他矩阵

      Q1 = R_cur * ego_size_ * R_cur.transpose();

      if (i == 0) {Q = Q1;}
      else
      {
        double beta = sqrt(Q1.trace() / Q2.trace());
        Q = (1 + 1 / beta) * Q1 + (1 + beta) * Q2;   /// Q:形状矩阵  安全椭球边界
      }


      Eigen::EigenSolver<Eigen::Matrix3d> es(Q);
      /// 计算 Q^(1/2）  real(): 实数类型
      Eigen::Matrix3d E = (es.eigenvectors() * (es.eigenvalues().cwiseSqrt()).asDiagonal() * es.eigenvectors().inverse()).real();
      /// 先检验是否满足安全上一个飞行通道约束，如果是，不需要重新为当前点生成走廊约束
      poly_indices(i) = getSikangConst(E); // push_back poly_constraints_
      ellipsoid_matrices_.push_back(E);

      /// 干扰椭球体
      Q2 = getDistrEllipsoid(Ts_, Q_init); //update Q_init  <1ms
    }
  }


  /// 更新力的结果
  void NMPCSolver::updateFORCESResults()
  {
    /* cmd update */
    pre_mpc_output_ = mpc_output_;
    have_mpc_traj_ = true;

    /// 这里只考虑了每个路径点的优化后航向在【-pi,pi】,参考的航向也是根据位置求的
    /// 和ego中航向处理稍微不一样（考虑最大航向转角）
    for (int i = 0; i < planning_horizon_; i++)
    {
      if (mpc_output_.at(i)(16) < -PI)
      {
        mpc_output_.at(i)(16) = mpc_output_.at(i)(16) + 2 * PI;
      }
      else if (mpc_output_.at(i)(16) > PI)
      {
        mpc_output_.at(i)(16) = mpc_output_.at(i)(16) - 2 * PI;
      }
    }

    mpc_output_.at(20) = mpc_output_.at(19);

    /* visulizations  */
    displayNMPCPoints();
    displayRefPoints();
    displayEllipsoids();
    displayPoly();

  }


  /// 根据欧拉角计算旋转矩阵    先旋转X轴再Y轴最后Z轴  从机体到惯性系
  Eigen::Matrix3d NMPCSolver::eulerToRot(Eigen::Vector3d &odom_euler)
  {

    Eigen::Vector3d rpy(odom_euler(0), odom_euler(1), odom_euler(2));
    Eigen::Quaternion<double> qx(cos(rpy(0) / 2), sin(rpy(0) / 2), 0, 0);
    Eigen::Quaternion<double> qy(cos(rpy(1) / 2), 0, sin(rpy(1) / 2), 0);
    Eigen::Quaternion<double> qz(cos(rpy(2) / 2), 0, 0, sin(rpy(2) / 2));
    Eigen::Matrix3d R = Eigen::Matrix3d(qz * qy * qx);

    return R;
  }

  /// 计算干扰椭球体
  // calculate the disturbance ellipsoid
  Eigen::Matrix3d NMPCSolver::getDistrEllipsoid(double t, Eigen::MatrixXd &Q_origin)
  {

    Eigen::MatrixXd Nt(nx, nx);
    Eigen::MatrixXd Array_Q(nx, nx);
    Eigen::MatrixXd temp_Q = Eigen::MatrixXd::Zero(nx, nx);
    double temp;
    Eigen::MatrixXd B = Phi_.transpose();   /// Phi_：系统矩阵
    //ROS_INFO_STREAM(" B : \n" << B );
    Eigen::ComplexSchur<Eigen::MatrixXd> SchurA(Phi_);  /// 舒尔分解 可以求非对称矩阵的特征值  A=UTU^T
    Eigen::MatrixXcd R = SchurA.matrixT();   /// 上三角矩阵
    Eigen::MatrixXcd U = SchurA.matrixU();   ///  酉矩阵
    Eigen::ComplexSchur<Eigen::MatrixXd> SchurB(B);
    Eigen::MatrixXcd S = SchurB.matrixT();
    Eigen::MatrixXcd V = SchurB.matrixU();
    // refer to this: https://stackoverflow.com/questions/56929966/implementing-the-bartels-stewart-algorithm-in-eigen3
    // ROS_INFO_STREAM(" R : \n" << R);
    // ROS_INFO_STREAM(" U : \n" << U );
    // ROS_INFO_STREAM(" S : \n" << S );
    // ROS_INFO_STREAM(" V : \n" << V );
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(nx, nx);

    for (int i = 0; i < nw; ++i)  /// nw = 3  nx = 状态变量数9
    {
      Nt = t * w_(i) * w_(i) * Dt_.col(i) * (Dt_.col(i)).transpose(); // nx * nx
      Array_Q = Nt - (-Phi_ * t).exp() * Nt * (-Phi_.transpose() * t).exp(); // nx * nx
      /// Phi*X + X*Phi' = W  求解李雅普诺夫方程：稠密的矩阵 Phi，Bartels-Stewart算法是最为高效的一种；对于大规模稀疏的矩阵 Phi ，若 W 是低秩矩阵，则Krylov类型的算法更加高效
      /// Bartels-Stewart算法
      Eigen::MatrixXcd F = (U.adjoint() * Array_Q) * V;  /// adjoint(): 伴随矩阵
      Eigen::MatrixXcd Y = Eigen::internal::matrix_function_solve_triangular_sylvester(R, S, F);
      Eigen::MatrixXd X = ((U * Y) * V.adjoint()).real();

      temp += sqrt(X.trace());
      temp_Q += X / sqrt(X.trace());
    }

    Eigen::MatrixXd Qd = temp * temp_Q; // nx * nx
    double beta = sqrt(Q_origin.trace() / Qd.trace());
    Eigen::MatrixXd Q_update = (1 + 1 / beta) * Q_origin + (1 + beta) * Qd; // nx * nx
    Eigen::MatrixXd position_Q = (Phi_ * t).exp() * Q_update * (Phi_.transpose() * t).exp();

    // update Q_origin
    Q_origin = Q_update;

    return position_Q.block(0, 0, 3, 3);
  }


  /// odom_euler:  roll pitch yaw
  /// 欧拉角
  Eigen::Matrix3d NMPCSolver::updateMatrix(Eigen::Vector3d &odom_euler, Eigen::Vector3d &odom_vel, double thrust_c)
  {
    // update At matrix
    double roll = odom_euler(0), pitch = odom_euler(1), yaw = odom_euler(2);
    double v1 = odom_vel(0), v2 = odom_vel(1), v3 = odom_vel(2);
    double comb0 = thrust_c * 1.0 / mass;

    double comb5 = cos(pitch)*sin(pitch);

    double comb6 = cos(pitch)*sin(roll);
    double comb7 = cos(pitch)*cos(roll);   ///
    double comb8 = sin(pitch)*cos(roll);
    double comb9 = sin(pitch)*sin(roll);
        
    double comb1 = cos(roll)*sin(yaw) - comb9*cos(yaw);
    double comb2 = sin(roll)*cos(yaw) - comb8*sin(yaw);    ///
    double comb3 = cos(roll)*cos(yaw) + comb9*sin(yaw);
    double comb4 = sin(roll)*sin(yaw) + comb8*cos(yaw);    ///

    // roll
    At_(3, 6) =  comb0*comb1;
    At_(4, 6) = -comb0*comb3;
    At_(5, 6) = -comb0*comb6;
    // pitch
    At_(3, 7) =  comb0*comb7*cos(yaw);
    At_(4, 7) =  comb0*comb7*sin(yaw);
    At_(5, 7) = -comb0*comb8;
    // yaw
    At_(3, 8) =  comb0*comb2;
    At_(4, 8) =  comb0*comb4;

    // add gradients with drag coefficient term (in some case we can neglect this term)
    Eigen::Matrix3d R_cur = eulerToRot(odom_euler);   ///

    At_.block(3, 3, 3, 3) = R_cur * drag_coefficient_matrix_ * R_cur.transpose(); //transpose does not change
    
    double drag = drag_coefficient_matrix_(0, 0);
    double cos_pitch_sq = cos(pitch)*cos(pitch);
    double sin_pitch_sq = sin(pitch)*sin(pitch);
    double cos_yaw_sq   = cos(yaw)*cos(yaw);
    double sin_yaw_sq   = sin(yaw)*sin(yaw);
    double sin_roll_sq  = sin(roll)*sin(roll);
    double temp_sqaure1 = pow(comb3, 2);
    double temp_sqaure2 = pow(comb1, 2);

    Eigen::Vector3d temp1, temp2;
    Eigen::Vector4d temp3;

    temp1(0) = comb6*comb4 - comb7*comb1;
    temp1(1) = comb3*comb4 + comb1*comb2;
    temp1(2) = comb6*comb2 - comb7*comb3;
 
    // roll
    At_(3, 6) +=  drag*(v3*temp1(0) + v2*temp1(1) - 2*v1*comb4*comb1);
    At_(4, 6) +=  drag*(v1*temp1(1) - v3*temp1(2) - 2*v2*comb3*comb2);
    At_(5, 6) +=  drag*(v1*temp1(0) - v2*temp1(2) + 2*v3*comb7*comb6);
 
    // pitch
    temp2(0) =   cos(yaw)*(sin_pitch_sq - cos_pitch_sq + cos_pitch_sq*sin_roll_sq) + comb9*comb1;
    temp2(1) =   2*comb5*cos(yaw)*sin(yaw) - comb6*(cos(yaw)*comb3 + sin(yaw)*comb1);
    temp2(2) =   sin(yaw)*(cos_pitch_sq - sin_pitch_sq - cos_pitch_sq*sin_roll_sq) + comb9*comb3;

    At_(3, 7) +=   drag*(v3*temp2(0) - v2*temp2(1) - v1*2*(comb5*cos_yaw_sq + comb6*comb1*cos(yaw)));
    At_(4, 7) += - drag*(v3*temp2(2) - v1*temp2(1) - v2*2*(comb5*sin_yaw_sq - comb6*comb3*sin(yaw)));
    At_(5, 7) +=   drag*(v1*temp2(0) - v2*temp2(2) + v3*2*(comb5 - comb5*sin_roll_sq));
 
    // yaw
    temp3(0) =   2*drag*(comb3*comb1 - cos_pitch_sq*cos(yaw)*sin(yaw));
    temp3(1) =   drag*(comb6*comb3 - comb5*sin(yaw));
    temp3(2) =   drag*(temp_sqaure1 - temp_sqaure2 - cos_pitch_sq*cos_yaw_sq + cos_pitch_sq*sin_yaw_sq);
    temp3(3) =   drag*(comb6*comb1 + comb5*cos(yaw));

    At_(3, 8) +=   v1*temp3(0)  - v3*temp3(1) - v2*temp3(2);
    At_(4, 8) += - v1*temp3(2)  - v3*temp3(3) - v2*temp3(0);
    At_(5, 8) += - v2*temp3(3)  - v1*temp3(1);
 
    // update Bt matrix
    Bt_(3, 3) =  1.0 / mass*comb4;
    Bt_(4, 3) = -1.0 / mass*comb2;
    Bt_(5, 3) =  1.0 / mass*comb7;

    Phi_ = At_ + Bt_ * Kt_;

    return R_cur;
  }

  /*visulizations*/
  void NMPCSolver::drawCmd(const Eigen::Vector3d &pos, const Eigen::Vector3d &vec, const int &id,
                           const Eigen::Vector4d &color)
  {
    visualization_msgs::Marker mk_state;
    mk_state.header.frame_id = "world";
    mk_state.header.stamp = ros::Time::now();
    mk_state.id = id;
    mk_state.type = visualization_msgs::Marker::ARROW;
    mk_state.action = visualization_msgs::Marker::ADD;

    mk_state.pose.orientation.w = 1.0;
    mk_state.scale.x = 0.1;
    mk_state.scale.y = 0.2;
    mk_state.scale.z = 0.3;

    geometry_msgs::Point pt;
    pt.x = pos(0);
    pt.y = pos(1);
    pt.z = pos(2);
    mk_state.points.push_back(pt);

    pt.x = pos(0) + vec(0);
    pt.y = pos(1) + vec(1);
    pt.z = pos(2) + vec(2);
    mk_state.points.push_back(pt);

    mk_state.color.r = color(0);
    mk_state.color.g = color(1);
    mk_state.color.b = color(2);
    mk_state.color.a = color(3);

    cmd_vis_pub_.publish(mk_state);
  }

  ///
  void NMPCSolver::displayNMPCPoints()
  {
    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.type = visualization_msgs::Marker::LINE_STRIP;
    mk.action = visualization_msgs::Marker::DELETE;

    nmpc_marker_pub_.publish(mk);  ///
    geometry_msgs::Point pt;
    std_msgs::ColorRGBA pc;

    mk.action = visualization_msgs::Marker::ADD;
    mk.pose.orientation.w = 1.0;
    mk.scale.x = 0.1;

    pc.r = 0.5;
    pc.g = 0;
    pc.b = 0;
    pc.a = 0.6;

    for (int i = 0; i < planning_horizon_; i++)
    {
      pt.x = pre_mpc_output_.at(i)(8);
      pt.y = pre_mpc_output_.at(i)(9);
      pt.z = pre_mpc_output_.at(i)(10);

      mk.points.push_back(pt);
      mk.colors.push_back(pc);
    }
    nmpc_marker_pub_.publish(mk);  ///
  }

  ///
  void NMPCSolver::displayRefPoints()
  {
    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.type = visualization_msgs::Marker::SPHERE_LIST;
    mk.action = visualization_msgs::Marker::DELETE;

    ref_marker_pub_.publish(mk); //
    geometry_msgs::Point pt;
    std_msgs::ColorRGBA pc;

    mk.action = visualization_msgs::Marker::ADD;
    mk.pose.orientation.w = 1.0;
    mk.scale.x = 0.2;
    mk.scale.y = 0.2;
    mk.scale.z = 0.2;

    pc.r = 0;
    pc.g = 0;
    pc.b = 0.5;
    pc.a = 0.6;

    for (int i = 0; i < planning_horizon_; i++)
    {
        /// kino搜索出的路径点
      pt.x = ref_total_pos_.at(i)(0);
      pt.y = ref_total_pos_.at(i)(1);
      pt.z = ref_total_pos_.at(i)(2);

      mk.points.push_back(pt);
      mk.colors.push_back(pc);
    }
    ref_marker_pub_.publish(mk); //
  }

  ///
  void NMPCSolver::displayEllipsoids()
  {
    decomp_ros_msgs::EllipsoidArray ellipsoids;
    for (int i = 0; i < planning_horizon_; i++)
    {
      decomp_ros_msgs::Ellipsoid ellipsoid;
      ellipsoid.d[0] = mpc_output_.at(i)(8);
      ellipsoid.d[1] = mpc_output_.at(i)(9);
      ellipsoid.d[2] = mpc_output_.at(i)(10);

      auto C = ellipsoid_matrices_[i];

      for (int x = 0; x < 3; x++)
      {
        for (int y = 0; y < 3; y++)
        {
          ellipsoid.E[3 * x + y] = C(x, y);
        }
      }
      ellipsoids.ellipsoids.push_back(ellipsoid);
    }

    ellipsoids.header.frame_id = "world";
    ellipsoid_pub_.publish(ellipsoids);  //
  }

  /// 根据前后位置计算偏航角
  void NMPCSolver::calculate_yaw(Eigen::Vector3d pos, Eigen::Vector3d pos_next)
  {
    double yaw = 0;

    Eigen::Vector3d dir = pos_next - pos;

    double yaw_temp = dir.norm() > 0.1 ? atan2(dir(1), dir(0)) : last_yaw_; // if dir is some or not change

    if (fabs(yaw_temp - last_yaw_) > PI)
    {
      if (yaw_temp > 0)
      { // -3.14 to 3.14
        yaw = yaw_temp - 2 * PI;
      }
      else
      { // 3.14 to -3.1
        yaw = yaw_temp + 2 * PI;
      }
    }
    else
    {
      yaw = yaw_temp;
    }

    yaw = 0.2 * last_yaw_ + 0.8 * yaw; // nieve LPF

    last_yaw_ = yaw;
    ref_yaw_ = yaw;
  }

  /*** callbacks  ***/
    void NMPCSolver::cmdTrajCallback(const ros::TimerEvent &e)
  {
    trajectory_msgs::MultiDOFJointTrajectory trajectory_msg;
    trajectory_msg.header.stamp = ros::Time::now();
    trajectory_msg.header.frame_id = "world";

    trajectory_msgs::MultiDOFJointTrajectoryPoint point_msg;
    point_msg.transforms.resize(1);
    point_msg.velocities.resize(1);
    point_msg.accelerations.resize(1);

    switch (cmd_status_)
    {
      case INIT_POSITION: /// if need any fixed position
      {
        break;
      }

      case ROTATE_YAW:  ///
      {
        /// init_yaw_是朝向目标点的方向，yaw_temp是根据的偏航角的角速度计算而来
        double yaw_temp = realOdom_(8) + (ros::Time::now() - change_yaw_time_).toSec() * init_yaw_dot_;
        double desired_yaw = init_yaw_ - realOdom_(8) >= 0 ? min(yaw_temp, init_yaw_) : max(yaw_temp, init_yaw_);
        /// 无人机悬停时旋转偏航角
        Eigen::Vector3d desired_position(realOdom_(0), realOdom_(1), realOdom_(2)); ///
        mav_msgs::msgMultiDofJointTrajectoryFromPositionYaw(desired_position, desired_yaw, &trajectory_msg);
        traj_cmd_pub_.publish(trajectory_msg);  ///
        break;
      }

      case WAIT:
      {
        break;
      }

      case PUB_TRAJ:
      {

        if (!have_mpc_traj_)  break;

        ros::Time time_now = ros::Time::now();
        /// pre_mpc_start_time_:没有考虑kino搜索时的时间
        double t_cur = (time_now - pre_mpc_start_time_).toSec();  /// 或者 t_cur = 当前时间 - 优化后的时间 - 设定的优化时间
        int cur_index = (int)(t_cur / Ts_);

        if (cur_index < planning_horizon_ - 1 && t_cur >= 0.0)
        {
            /// 线性插值
          Eigen::VectorXd mpcq = pre_mpc_output_.at(cur_index) + fmod(t_cur, Ts_) / Ts_ * (pre_mpc_output_.at(cur_index + 1) - pre_mpc_output_.at(cur_index));

          /// 期望的位置
          point_msg.transforms[0].translation.x = mpcq(8);
          point_msg.transforms[0].translation.y = mpcq(9);
          point_msg.transforms[0].translation.z = mpcq(10);

          /// 期望的速度
          point_msg.velocities[0].linear.x = mpcq(11);
          point_msg.velocities[0].linear.y = mpcq(12);
          point_msg.velocities[0].linear.z = mpcq(13);

          /// 都是0
          point_msg.velocities[0].angular.x = mpcq(0);
          point_msg.velocities[0].angular.y = mpcq(1);
          point_msg.velocities[0].angular.z = mpcq(2);


          Eigen::Vector3d odom_euler(mpcq(14), mpcq(15), mpcq(16));
          Eigen::Vector3d thrust_b(0.0, 0.0, mpcq(3));
          Eigen::Vector3d thrust_w = eulerToRot(odom_euler) * thrust_b;

          /// 期望的就速度
          point_msg.accelerations[0].linear.x = thrust_w(0) / mass;
          point_msg.accelerations[0].linear.y = thrust_w(1) / mass;
          point_msg.accelerations[0].linear.z = thrust_w(2) / mass - g;

          /// 四元数 ：按XYZ的顺序左乘的
          Eigen::Quaterniond q = Eigen::AngleAxisd(mpcq(14), Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(mpcq(15), Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(mpcq(16), Eigen::Vector3d::UnitZ());

          point_msg.transforms[0].rotation.x = q.x();
          point_msg.transforms[0].rotation.y = q.y();
          point_msg.transforms[0].rotation.z = q.z();
          point_msg.transforms[0].rotation.w = q.w();

          Eigen::Vector3d pos(mpcq(8), mpcq(9), mpcq(10));
          Eigen::Vector3d dir(cos(mpcq(16)), sin(mpcq(16)), 0.0);
          drawCmd(pos, dir, 2, Eigen::Vector4d(1, 1, 0, 0.7));

          trajectory_msg.points.push_back(point_msg);
          traj_cmd_pub_.publish(trajectory_msg);
          finish_mpc_cmd_ = false;
        }
        else
        {
          finish_mpc_cmd_ = true;
          if (pub_end_)
              cmd_status_ = CMD_STATUS::PUB_END;
        }
        break;
      }

      case PUB_END:
      {
        point_msg.transforms[0].translation.x = cmd_end_pt_(0);
        point_msg.transforms[0].translation.y = cmd_end_pt_(1);
        point_msg.transforms[0].translation.z = cmd_end_pt_(2);

        point_msg.velocities[0].linear.x = 0;
        point_msg.velocities[0].linear.y = 0;
        point_msg.velocities[0].linear.z = 0;

        double roll  = pre_mpc_output_.at(19)(14);
        double pitch = pre_mpc_output_.at(19)(15);
        double yaw   = pre_mpc_output_.at(19)(16);

        ///
        Eigen::Quaterniond q = Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ());

        point_msg.transforms[0].rotation.x = q.x();
        point_msg.transforms[0].rotation.y = q.y();
        point_msg.transforms[0].rotation.z = q.z();
        point_msg.transforms[0].rotation.w = q.w();

        trajectory_msg.points.push_back(point_msg);

        traj_cmd_pub_.publish(trajectory_msg);
        initialized_output_ = false;    /// 初始化没有完成
        ROS_INFO_STREAM(" end time is " << ros::Time::now().toSec());
        cmd_status_ = CMD_STATUS::WAIT;  /// 到达终点 切换等待
        break;
      }
    }
  }

  //pcl::PointCloud< pcl::PointXYZ > pcl的点云
  //pcl::PCLPointCloud2 pcl的第⼆种点云
  //sensor_msgs::PointCloud2 ROS中的点云
  //sensor_msgs::PointCloud ROS中的点云
  /// 有的传感器驱动数据格式就是PointCloud2
  void NMPCSolver::cloudCallback(const sensor_msgs::PointCloud2 &msg)
  {
    sensor_msgs::PointCloud out_cloud;
    sensor_msgs::convertPointCloud2ToPointCloud(msg, out_cloud);  /// 先进行转换PointCloud
    vec_obs_ = DecompROS::cloud_to_vec(out_cloud);  /// vec_obs_:点云即是障碍物
  }

}
