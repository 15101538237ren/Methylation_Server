def fake_nearby_simulation(function_util):
    random_propensity_list = function_util.set_collaborative_params(U_plus_in=0.005,H_plus_in=0.004,M_minus_in=0.037,H_minus_in=0.034,H_p_H_in=0.233,H_p_M_in=0.233,U_p_M_in=0.232,H_m_U_in=0.084,M_m_U_in=0.084)

    sim_rounds = range(1,3)

    nearby_index = "fake_nearby_simulation_try2"
    NEARBY_OUTPUT_DIR = "data" + os.sep + nearby_index

    max_cpg_sites = 200
    generations = 10
    multi_threads = True
    nearby = 1
    real_nearby=False
    n_time_step=100


    max_cells = 2
    detail_for_timestep = range(99)

    geometric_p = 0.3
    plot = False
    cpg_max_pos, pos_list = function_util.construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p,
                                                                                     plot=plot)
    m_ratio = 0.181214  # the site origin ratio
    h_ratio = 0.427782
    u_ratio = 0.391004

    init_cell = function_util.generate_CpG_in_methylation_percent_UHM(max_cpg_sites, m_ratio,u_ratio)  # generate a cpg chain which have the methylation status
    simulator = Simulator(random_propensity_list, rounds=sim_rounds, out_dir=NEARBY_OUTPUT_DIR,
                      max_cpg_sites=max_cpg_sites, generations=generations, pos_list=pos_list,
                      multi_threads=multi_threads, init_cell=init_cell, nearby=nearby, max_cells=max_cells,
                      index=nearby_index, detail_for_timestep=detail_for_timestep,real_nearby=real_nearby,n_time_step=n_time_step)
    simulator.run()

    sorted_ratio_dir_name = "sorted_ratio"
    sorted_ratio_bk_dir_name = "sorted_ratio_bk"
    sorted_detail_dir_name = "sorted_detail"
    bed_files_dir_name = "bed_files"
    rd_with_dir_name = "rd_with"
    rd_without_dir_name = "rd_without"

    sorted_ratio_dir = NEARBY_OUTPUT_DIR + os.sep + sorted_ratio_dir_name
    sorted_ratio_bk_dir = NEARBY_OUTPUT_DIR + os.sep + sorted_ratio_bk_dir_name
    sort_detail_dir = NEARBY_OUTPUT_DIR + os.sep + sorted_detail_dir_name
    bed_files_dir = NEARBY_OUTPUT_DIR + os.sep + bed_files_dir_name
    rd_with_dir = NEARBY_OUTPUT_DIR + os.sep + rd_with_dir_name
    rd_without_dir = NEARBY_OUTPUT_DIR + os.sep + rd_without_dir_name

    simulator.sort_the_simulaiton_result(NEARBY_OUTPUT_DIR,sorted_ratio_dir,sort_detail_dir, simulator.rounds[0],simulator.rounds[len(simulator.rounds)-1], [], False)
    shutil.copytree(sorted_ratio_dir, sorted_ratio_bk_dir)
    filter_bounds = simulator.set_filter_range_bounds(m_down=0.10, m_up=0.28, h_down=0.30, h_up=0.55, u_down=0.30,
                                                      u_up=0.55)

    remained_generations = simulator.sort_the_simulaiton_result(NEARBY_OUTPUT_DIR, sorted_ratio_dir,
                                                                sort_detail_dir, simulator.rounds[0],
                                                                simulator.rounds[len(simulator.rounds) - 1], [], True,
                                                                filter_bounds)  # filter the sort result according to the filter bound for m,h,u ratio
    simulator.sort_to_bed(simulator.pos_list,sort_detail_dir,bed_files_dir,gens=remained_generations)

    calc_d_max = 50  # 计算的相关性最大距离
    calc_interval = False  # 是否包含中间的位点
    ignore_d = False  # 是否忽略位点间距离而计算相关性
    list_gen = remained_generations
    str_list_gen = []
    for item in list_gen:
        item = str(item).replace(".", "_")
        str_list_gen.append(str(item))

    simulator.calc_corr(bed_files_dir, rd_with_dir, rd_without_dir, str_list_gen, calc_d_max,
                        calc_interval=calc_interval, ignore_d=ignore_d)
def random_col_simulation(function_util):
    random_propensity_list = function_util.set_collaborative_params(U_plus_in=0.0005,H_plus_in=0.0005,M_minus_in=0.07,H_minus_in=0.069 ,H_p_H_in=0.236,H_p_M_in=0.236,U_p_M_in=0.236 ,H_m_U_in=0.06,M_m_U_in=0.06)

    sim_rounds = range(1,3)

    random_index = "random_simulation_try4"
    RANDOM_OUTPUT_DIR = "data" + os.sep + random_index

    max_cpg_sites = 200
    generations =1000
    multi_threads = True
    nearby = -1
    max_cells = 2
    detail_for_timestep = [0, 1, 2]
    excepts_round=[]
    geometric_p = 0.3
    plot = False
    cpg_max_pos, pos_list = function_util.construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p,
                                                                                     plot=plot)
    m_ratio = 0.181214  # the site origin ratio
    h_ratio = 0.427782
    u_ratio = 0.391004

    init_cell = function_util.generate_CpG_in_methylation_percent_UHM(max_cpg_sites, m_ratio,u_ratio)  # generate a cpg chain which have the methylation status
    simulator = Simulator(random_propensity_list, rounds=sim_rounds, out_dir=RANDOM_OUTPUT_DIR,
                      max_cpg_sites=max_cpg_sites, generations=generations, pos_list=pos_list,
                      multi_threads=multi_threads, init_cell=init_cell, nearby=nearby, max_cells=max_cells,
                      index=random_index, detail_for_timestep=detail_for_timestep)
    simulator.run()

    sorted_ratio_dir_name = "sorted_ratio"
    sorted_ratio_bk_dir_name = "sorted_ratio_bk"
    sorted_detail_dir_name = "sorted_detail"
    bed_files_dir_name = "bed_files"
    rd_with_dir_name = "rd_with"
    rd_without_dir_name = "rd_without"

    sorted_ratio_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_ratio_dir_name
    sorted_ratio_bk_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_ratio_bk_dir_name
    sort_detail_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_detail_dir_name
    bed_files_dir = RANDOM_OUTPUT_DIR + os.sep + bed_files_dir_name
    rd_with_dir = RANDOM_OUTPUT_DIR + os.sep + rd_with_dir_name
    rd_without_dir = RANDOM_OUTPUT_DIR + os.sep + rd_without_dir_name

    simulator.sort_the_simulaiton_result(RANDOM_OUTPUT_DIR,sorted_ratio_dir,sort_detail_dir, simulator.rounds[0],simulator.rounds[len(simulator.rounds)-1],excepts_round, False)
    shutil.copytree(sorted_ratio_dir, sorted_ratio_bk_dir)
    filter_bounds = simulator.set_filter_range_bounds(m_down=0.10, m_up=0.28, h_down=0.30, h_up=0.55, u_down=0.30,
                                                      u_up=0.55)

    remained_generations = simulator.sort_the_simulaiton_result(RANDOM_OUTPUT_DIR, sorted_ratio_dir,
                                                                sort_detail_dir, simulator.rounds[0],
                                                                simulator.rounds[len(simulator.rounds) - 1], excepts_round, True,
                                                                filter_bounds)  # filter the sort result according to the filter bound for m,h,u ratio
    simulator.sort_to_bed(simulator.pos_list,sort_detail_dir,bed_files_dir,gens=remained_generations)

    calc_d_max = 50  # 计算的相关性最大距离
    calc_interval = False  # 是否包含中间的位点
    ignore_d = False  # 是否忽略位点间距离而计算相关性
    list_gen = remained_generations
    str_list_gen = []
    for item in list_gen:
        item = str(item).replace(".", "_")
        str_list_gen.append(str(item))

    simulator.calc_corr(bed_files_dir, rd_with_dir, rd_without_dir, str_list_gen, calc_d_max,
                        calc_interval=calc_interval, ignore_d=ignore_d)
def traditional_simulation(function_util):
    # traditional simulation
    traditional_propensity_list = function_util.set_standard_params(U_plus_in=0.05, H_plus_in=0.05, M_minus_in=0.05,
                                                                    H_minus_in=0.05)
    sim_rounds = range(1,3)

    traditional_index = "traditional_simulation"
    TRADITIONAL_OUTPUT_DIR = "data" + os.sep + traditional_index

    max_cpg_sites = 500
    generations = 10
    multi_threads = True
    nearby = -1
    max_cells = 2
    detail_for_timestep = [0, 1, 2]

    geometric_p = 0.3
    plot = False
    cpg_max_pos, pos_list = function_util.construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p,
                                                                                     plot=plot)
    m_ratio = 0.181214  # the site origin ratio
    h_ratio = 0.427782
    u_ratio = 0.391004

    init_cell = function_util.generate_CpG_in_methylation_percent_UHM(max_cpg_sites, m_ratio,
                                                                      u_ratio)  # generate a cpg chain which have the methylation status
    simulator = Simulator(traditional_propensity_list, rounds=sim_rounds, out_dir=TRADITIONAL_OUTPUT_DIR,
                          max_cpg_sites=max_cpg_sites, generations=generations, pos_list=pos_list,
                          multi_threads=multi_threads, init_cell=init_cell, nearby=nearby, max_cells=max_cells,
                          index=traditional_index, detail_for_timestep=detail_for_timestep)
    simulator.run()

    sorted_ratio_dir_name="sorted_ratio"
    sorted_ratio_bk_dir_name="sorted_ratio_bk"
    sorted_detail_dir_name="sorted_detail"
    bed_files_dir_name="bed_files"
    rd_with_dir_name="rd_with"
    rd_without_dir_name = "rd_without"

    sorted_ratio_dir = TRADITIONAL_OUTPUT_DIR + os.sep + sorted_ratio_dir_name
    sorted_ratio_bk_dir = TRADITIONAL_OUTPUT_DIR + os.sep + sorted_ratio_bk_dir_name
    sort_detail_dir = TRADITIONAL_OUTPUT_DIR + os.sep + sorted_detail_dir_name
    bed_files_dir = TRADITIONAL_OUTPUT_DIR + os.sep + bed_files_dir_name
    rd_with_dir = TRADITIONAL_OUTPUT_DIR + os.sep + rd_with_dir_name
    rd_without_dir = TRADITIONAL_OUTPUT_DIR + os.sep + rd_without_dir_name

    simulator.sort_the_simulaiton_result(TRADITIONAL_OUTPUT_DIR,sorted_ratio_dir,sort_detail_dir, simulator.rounds[0],simulator.rounds[len(simulator.rounds)-1], [], False)
    shutil.copytree(sorted_ratio_dir, sorted_ratio_bk_dir)
    filter_bounds=simulator.set_filter_range_bounds(m_down=0.10,m_up=0.24,h_down=0.30,h_up=0.52,u_down=0.30,u_up=0.52)

    remained_generations=simulator.sort_the_simulaiton_result(TRADITIONAL_OUTPUT_DIR,sorted_ratio_dir,sort_detail_dir, simulator.rounds[0],simulator.rounds[len(simulator.rounds)-1], [], True,filter_bounds) # filter the sort result according to the filter bound for m,h,u ratio
    simulator.sort_to_bed(simulator.pos_list,sort_detail_dir,bed_files_dir,gens=remained_generations)

    calc_d_max=500 #计算的相关性最大距离
    calc_interval=False #是否包含中间的位点
    ignore_d=False #是否忽略位点间距离而计算相关性
    list_gen = remained_generations
    str_list_gen = []
    for item in list_gen:
        item = str(item).replace(".", "_")
        str_list_gen.append(str(item))

    simulator.calc_corr(bed_files_dir,rd_with_dir,rd_without_dir,str_list_gen,calc_d_max,calc_interval=calc_interval,ignore_d=ignore_d)
