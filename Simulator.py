# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections
import numpy as np
from FunctionUtil import FunctionUtil

BED_FILE_PATH="chr1.bed"
THRESHOLD_U=0.2
THRESHOLD_M=0.8
N_STEP=100

REACTION_STATUS_HASH={0:'H',1:'M',2:'H',3:'U',4:'M',5:'M',6:'H',7:'U',8:'H'}
RIGHT_STATUS_OF_REACTION={0:"U",1:"H",2:"M",3:"H",4:"H",5:"H",6:"U",7:"H",8:"M"}
STATE_OF_COLLABOR_REACTION={4:"H",5:"M",6:"M",7:"U",8:"U"}
BASE_RATE_HASH={4:1,5:1,6:0,7:3,8:2} #collaboration reaction base reaction rate index_hash
COLLABORATION_INDEX=4
class Simulator(object):
    '''
        This is a base class for nearby , random collaborative and traditional simulation
    '''

    # Construction function
    def __init__(self,propensity_list,rounds=range(1,2),out_dir="out",max_cpg_sites=1000,generations=10,pos_list=[],multi_threads=False,init_cell="",nearby=-1,max_cells=2,index="index",detail_for_timestep=[0,1]):
        #create the output dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.out_dir=out_dir
        self.propensity_list=propensity_list
        self.max_cpg_sites=max_cpg_sites
        self.generations=generations
        self.rounds=rounds
        self.pos_list=pos_list
        self.multi_threads=multi_threads
        self.init_cell=init_cell
        self.nearby=nearby
        self.max_cells=max_cells
        self.index=index
        self.lock = threading.Lock()
        self.detail_for_timestep=detail_for_timestep
        if self.multi_threads==True:
            self.num_sites_per_thread=100
            self.max_threads = 250
            self.num_of_threads = int(math.ceil(float(len(init_cell)) / self.num_sites_per_thread)) # thread counts
            self.num_turns = int(math.ceil(float(self.num_of_threads) / self.max_threads)) # when the max_thread is not enough for one round simulation, split it in to num_turns
            self.threads = [] # thread array
            print "%d thread is needed!" % self.num_of_threads
    def run(self):
        starttime_new=datetime.datetime.now()
        for i_round in self.rounds:
            self.threads=[]
            starttime_a_round = datetime.datetime.now()

            detail_file_full_path = self.out_dir + os.sep + "detail_" + str(i_round) + ".csv" # the sites status detail statistics file
            detail_file = open(detail_file_full_path, "w")

            ratio_file_full_path = self.out_dir + os.sep + "ratio_" + str(i_round) + ".csv" # the sites status ratio statistics file
            ratio_file = open(ratio_file_full_path, "w")

            if self.multi_threads == True: #if multi threads
                for thread_no in range(self.num_of_threads): # create all threads
                    start = thread_no * self.num_sites_per_thread
                    if thread_no != self.num_of_threads - 1: # get the init cell status
                        end=(thread_no + 1) * self.num_sites_per_thread
                    else:
                        end=len(self.init_cell)

                    init_cell_of_thread = self.init_cell[start:end]
                    pos_tmp_list = self.pos_list[start:end]
                    # create thread
                    th = threading.Thread(target=self.multi_thread_simulation, args=(
                    i_round, thread_no, self.generations, N_STEP, init_cell_of_thread, detail_file,self.detail_for_timestep , ratio_file,
                    self.propensity_list, self.index,self.nearby, self.max_cells, pos_tmp_list))
                    self.threads.append(th)

                for i in range(self.num_turns):
                    if i != self.num_turns - 1:
                        threads_to_simulate = self.threads[i * self.max_threads:(i + 1) * self.max_threads]
                    else:
                        threads_to_simulate = self.threads[i * self.max_threads:len(self.threads)]
                    # 启动线程
                    for t in threads_to_simulate:
                        try:
                            t.start()
                        except RuntimeError,e:
                            print e
                    for t in threads_to_simulate:
                        t.join()
                    print '%d round:%d turn thread is waitting for next execution...' % (i_round, i)
                detail_file.close()
                ratio_file.close()
            else:
                self.multi_thread_simulation(i_round,0,self.generations, N_STEP, self.init_cell,
                                  detail_file, self.detail_for_timestep, ratio_file, self.propensity_list, self.index, nearby=self.nearby, max_cells= self.max_cells,index_pos_list=self.pos_list)
            endtime_a_round = datetime.datetime.now()
            print "One Round in " + str((endtime_a_round - starttime_a_round).seconds) + " seconds\n"
        endtime = datetime.datetime.now()
        print "running "+str(self.rounds)+" rounds simulation in " + str((endtime - starttime_new).seconds) + " s\n"
    def multi_thread_write_of_a_list(self,file,list):
        self.lock.acquire()
        for item in list:
            print >> file, item
        self.lock.release()
    def multi_thread_simulation(self, times_idx, thread_no, generations, n_time_step, init_cell, detail_file, detail_for_time_steps, ratio_file, propensity_list, exp_name, nearby=-1, max_cells=1000, index_pos_list=[]):
        print "%d round: Thread %d is started!" % (times_idx, thread_no)
        cell_collection = [init_cell]

        for i in range(1, generations + 1): # loop the generations
            cells_wait_to_add = []

            if len(cell_collection) > max_cells / 2:
                cell_collection = random.sample(cell_collection, max_cells / 2)

            time_old = datetime.datetime.now()
            print "%d round: thread %d : %s  gen: %d, father col len:%d, sons col len:%d now started!" % (
            times_idx, thread_no, exp_name, i, len(cell_collection), len(cells_wait_to_add))

            M_count_statistics = []
            H_count_statistics = []
            U_count_statistics = []

            for j in range(n_time_step):
                M_count_statistics.append([])
                H_count_statistics.append([])
                U_count_statistics.append([])
            M_count_statistics, H_count_statistics, U_count_statistics, out_detail_seq_arr, cell_collection, cells_wait_to_add = self.simulate_common(
                cell_collection, n_time_step, propensity_list, M_count_statistics, H_count_statistics,
                U_count_statistics, cells_wait_to_add, nearby, detail_for_time_steps,
                index_pos_list=index_pos_list)

            m_means_ratio = (np.mean(np.array(M_count_statistics), axis=1) / len(init_cell)).tolist()
            h_means_ratio = (np.mean(np.array(H_count_statistics), axis=1) / len(init_cell)).tolist()
            u_means_ratio = (np.mean(np.array(U_count_statistics), axis=1) / len(init_cell)).tolist()

            print "%d round: thread %d finished and starting to write statistics file!" % (times_idx, thread_no)
            out_str_of_statistics = []
            for t in range(n_time_step):
                if t in detail_for_time_steps:
                    time = (i - 1) + t / float(n_time_step)
                    out_str = str(thread_no) + "," + str(round(time, 2)) + "," + str(
                        round(m_means_ratio[t], 4)) + "," + str(round(h_means_ratio[t], 4)) + "," + str(
                        round(u_means_ratio[t], 4)) + "\n"
                    out_str_of_statistics.append(out_str)

            self.multi_thread_write_of_a_list(ratio_file, out_str_of_statistics) # write the ratio file in the multi_thread env
            print "%d round:thread %d finished writing statistics file and start to wrt detail file!" % (
            times_idx, thread_no)

            detail_seq_str_arr = []
            for idx, item in enumerate(out_detail_seq_arr):
                for j, detail_seq in enumerate(item):
                    time = (i - 1) + detail_for_time_steps[j] / float(n_time_step)
                    out_str = str(thread_no) + "," + str(round(time, 2)) + "," + str(idx) + "," + detail_seq
                    detail_seq_str_arr.append(out_str)

            self.multi_thread_write_of_a_list(detail_file, detail_seq_str_arr)
            print "%d round:thread %d finished writing detail file!" % (times_idx, thread_no)

            time_collapsed = str((datetime.datetime.now() - time_old).seconds)
            print "%d round:thread %d : %s  gen: %d, father col len:%d, sons col len:%d in %s seconds" % (
            times_idx, thread_no, exp_name, i, len(cell_collection), len(cells_wait_to_add), time_collapsed)

            cell_collection = cells_wait_to_add
    def simulate_common(self,cell_collection,n_time_step,PROPENCITY_LIST,M_count_statistics,H_count_statistics,U_count_statistics,cells_wait_to_add,nearby,detail_for_time_steps=[],index_pos_list=[]):
        out_detail_seq_arr=[]
        for idx,cell in enumerate(cell_collection): #loop cell in cell_collection

            out_detail_seq_arr.append([])

            for j in range(n_time_step): #loop the time step

                for k in range(len(cell)): #loop for every site in a cell
                    #get the reaction site index:target_reaction_CpG_site ,and collaborative site: col_CpG_site_index

                    target_reaction_CpG_site=random.randint(0,len(cell_collection[idx])-1)

                    if nearby>0:
                        start=max(target_reaction_CpG_site-nearby,0)
                        end=min(target_reaction_CpG_site+nearby,len(cell_collection[idx])-1)
                    else:
                        start=0
                        end=len(cell_collection[idx])-1

                    col_CpG_site_index=random.randint(start,end)
                    while(target_reaction_CpG_site==col_CpG_site_index):
                        col_CpG_site_index=random.randint(start,end)

                    #最后根据距离算propencity_list中对应的反应概率
                    if len(index_pos_list)>0:
                        # 算出两者对应于index的距离
                        pos_target = index_pos_list[target_reaction_CpG_site]
                        col_site_pos = index_pos_list[col_CpG_site_index]

                        distance=int(math.fabs(pos_target-col_site_pos))
                        reaction_ratio=self.phi(d=distance)
                        propensity_tmp=self.scale_propensity_list(reaction_ratio,PROPENCITY_LIST)
                    else:
                        propensity_tmp=PROPENCITY_LIST

                    #print_propensity_list(propensity_tmp,collaborative=True)
                    sum_propensity = 0.0
                    for item in propensity_tmp:
                        sum_propensity = sum_propensity + item

                    num_of_reactions = len(propensity_tmp)
                    random_number = random.random()

                    #根据反应速率由Gillespie选择一个反应
                    reaction_id = self.select_reaction(propensity_tmp, num_of_reactions, sum_propensity, random_number)

                    #最后完成反应后该细胞的该位点的状态
                    status_of_target_site=cell_collection[idx][target_reaction_CpG_site]

                    #目标位点能否发生这个反应，不能则直接continue当前位点的尝试
                    if RIGHT_STATUS_OF_REACTION[reaction_id]!=status_of_target_site:
                        continue
                    #若当前位点是协作反应
                    if reaction_id >= COLLABORATION_INDEX:
                        #协作反应的状态
                        status_of_col_site=cell_collection[idx][col_CpG_site_index]
                        #若协作位点满足协作反应的位点条件
                        if STATE_OF_COLLABOR_REACTION[reaction_id]==status_of_col_site:
                            CpG_str_list=list(cell_collection[idx])
                            CpG_str_list[target_reaction_CpG_site]=REACTION_STATUS_HASH[reaction_id]
                            cell_collection[idx]=''.join(CpG_str_list)
                        else:
                            continue
                    else:
                        #非协作反应,孤立反应
                        cell_status=cell_collection[idx]
                        CpG_pre_str=cell_status[0:target_reaction_CpG_site]
                        CpG_post_str = cell_status[target_reaction_CpG_site+1:len(cell_status)]
                        cell_collection[idx]=CpG_pre_str+REACTION_STATUS_HASH[reaction_id]+CpG_post_str
                        #print "pre is %s and post is %s" %(cell_status[target_reaction_CpG_site],cell_collection[idx][target_reaction_CpG_site])
                if M_count_statistics!=None and H_count_statistics!=None and U_count_statistics!=None:
                    m_count=cell_collection[idx].count("M")
                    M_count_statistics[j].append(m_count)
                    h_count=cell_collection[idx].count("H")
                    H_count_statistics[j].append(h_count)
                    u_count=cell_collection[idx].count("U")
                    U_count_statistics[j].append(u_count)
                if len(detail_for_time_steps)!=0:
                    if j in detail_for_time_steps:
                        out_detail_seq_arr[idx].append(cell_collection[idx])
            cell_of_source=cell_collection[idx]
            [cell1,cell2]=self.cell_division(cell_of_source)
            cells_wait_to_add.append(cell1)
            cells_wait_to_add.append(cell2)
        if M_count_statistics!=None and H_count_statistics!=None and U_count_statistics!=None:
            return M_count_statistics,H_count_statistics,U_count_statistics,out_detail_seq_arr,cell_collection,cells_wait_to_add
        else:
            return cell_collection,cells_wait_to_add
    def phi(self,d=2): #the phi function which used for control the collaborative rate
        return 1.0
    def scale_propensity_list(self,ratio,propensity_list): # according to the ratio to scale the propensity_list collaborative rate
        if len(propensity_list) <= 4:
            return propensity_list
        else:
            propensity_list_new = []
            for i in range(0,9):
                if i >= 4:
                    base_rate = propensity_list[BASE_RATE_HASH[i]] # calculate the reaction rate
                    scale = propensity_list[i] - base_rate
                    new_reaction_rate = base_rate + scale * ratio

                    propensity_list_new.append(new_reaction_rate)
                else:
                    propensity_list_new.append(propensity_list[i])  # add in the non-collaborative reaction rate
            return propensity_list_new
    def select_reaction(self,propencity_list, num_of_reactions, sum_propencity, random_number):
        reaction = -1
        tmp_sum_propencity = 0.0
        random_number = random_number * sum_propencity
        for i in range(num_of_reactions + 1):
            tmp_sum_propencity = tmp_sum_propencity + propencity_list[i]
            if random_number < tmp_sum_propencity:
                reaction = i
                break
        return reaction
    def cell_division(self,cell):
        cell1 = ""
        cell2 = ""
        for idx, val in enumerate(cell):
            if val == "M":
                cell1 = cell1 + "H"
                cell2 = cell2 + "H"
            elif val == "U":
                cell1 = cell1 + "U"
                cell2 = cell2 + "U"
            elif val == "H":
                if random.random() > 0.5:
                    cell1 = cell1 + "H"
                    cell2 = cell2 + "U"
                else:
                    cell1 = cell1 + "U"
                    cell2 = cell2 + "H"
        return [cell1, cell2]
    
    def sort_the_simulaiton_result(self,base_dir,sorted_ratio_dir,sort_detail_dir,start_round,end_round,excepts_rounds=[],is_filter=False,filter_bounds=[]):# 对多线程的模拟结果进行排序整合成按照generation的统计文件,见sorted_ratio/和sorted_detail文件夹
        if not os.path.exists(sorted_ratio_dir):
            os.makedirs(sorted_ratio_dir)

        if not os.path.exists(sort_detail_dir):
            os.makedirs(sort_detail_dir)
        remained_generations=self.sort_ratio(base_dir,sorted_ratio_dir,start_round,end_round,excepts_rounds=excepts_rounds,is_filter=is_filter,filter_bounds=filter_bounds)
        self.sort_detail(base_dir,sort_detail_dir,start_round,end_round,excepts_rounds)

        return remained_generations

    def set_filter_range_bounds(self,m_down=0.10,m_up=0.24,h_down=0.30,h_up=0.52,u_down=0.30,u_up=0.52):
        return [m_down,m_up,h_down,h_up,u_down,u_up]
    def sort_ratio(self,input_file_dir, output_file_dir, start_round, end_round, excepts_rounds=[], is_filter=False,filter_bounds=[]): # 对各个线程产生的ratio进行排序,生成gen_X_XX_ratio.csv
        pattern = r'(\d+),([\d]+\.[\d]*),([\d]+\.[\d]*),([\d]+\.[\d]*),([\d]+\.[\d]*)'
        hash_of_sort_ratio = {}
        for i in range(start_round, end_round + 1):
            if i in excepts_rounds:
                continue
            hash_of_i = {}
            infile_path = input_file_dir + os.sep + "ratio_" + str(i) + ".csv"
            infile = open(infile_path, "r")
            line = infile.readline()
            while line:
                match = re.search(pattern, line)
                if match:
                    thread_no = int(match.group(1))
                    generation_step = float(match.group(2))
                    m_ratio = float(match.group(3))
                    h_ratio = float(match.group(4))
                    u_ratio = float(match.group(5))

                    ratio_list = [m_ratio, h_ratio, u_ratio]

                    if generation_step not in hash_of_i.keys():
                        hash_of_i[generation_step] = {}
                    hash_of_i[generation_step][thread_no] = ratio_list

                line = infile.readline()
                if not line:
                    break
            hash_of_i_new = {}
            for (key, val) in hash_of_i.items():
                sorted_hash_of_gen = sorted(val.iteritems(), key=lambda d: d[0])
                hash_of_i_new[key] = sorted_hash_of_gen
            hash_of_sort_ratio[i] = hash_of_i_new
            infile.close()
        generation_steps = hash_of_sort_ratio[start_round].keys()
        remained_generations=[]
        for gen_step in generation_steps:
            gen_step_str = str(gen_step).replace(".", "_")
            out_file_path = output_file_dir + os.sep + "gen_" + gen_step_str + "_ratio.csv"
            out_file = open(out_file_path, 'w')
            gen_is_good = False
            for round_i in range(start_round, end_round + 1):
                if round_i in excepts_rounds:
                    continue
                seq_list = hash_of_sort_ratio[round_i][gen_step]

                # 1.产生比例的时候把此部分取消注释
                if is_filter == False:
                    out_file.write(str(round_i) + ",")

                    len_item_of_seq_list = len(seq_list[0][1])
                    thread_num = len(seq_list)
                    for m_h_u_i in range(len_item_of_seq_list):
                        for thread_no in range(thread_num):
                            out_file.write(str(seq_list[thread_no][1][m_h_u_i]) + ",")
                    out_file.write("\n")
                # 1.end

                # 2.统计可用的参数时把此部分取消注释
                elif is_filter == True:
                    m_sum = 0
                    u_sum = 0
                    h_sum = 0

                    index = 0
                    thread_num = len(seq_list)
                    for thread_no in range(thread_num):
                        index = index + 1
                        m_sum = m_sum + seq_list[thread_no][1][0]
                        h_sum = h_sum + seq_list[thread_no][1][1]
                        u_sum = u_sum + seq_list[thread_no][1][2]
                    m_mean = m_sum / float(index)
                    h_mean = h_sum / float(index)
                    u_mean = u_sum / float(index)

                    if m_mean > filter_bounds[0] and m_mean < filter_bounds[1] and h_mean > filter_bounds[2] and h_mean < filter_bounds[3] and u_mean > filter_bounds[4] and u_mean < filter_bounds[5]:
                        print "gen: %f round:%d ,m:%f ,h:%f ,u:%f" % (gen_step, round_i, m_mean, h_mean, u_mean)
                        gen_is_good = True
                        out_file.write(str(round_i) + "," + str(m_mean) + "," + str(h_mean) + "," + str(u_mean) + "\n")
                        if gen_step not in remained_generations:
                            remained_generations.append(gen_step)
                        # 2.end
            out_file.close()
            if gen_is_good==False and is_filter==True:
                os.remove(out_file_path)
        if is_filter==True:
            return remained_generations
    def sort_detail(self,input_file_dir, output_file_dir,start_round,end_round, excepts_rounds=[]): # 对各个线程产生的位点状态结果进行排序,生成一个文件gen_X_XX_detail.csv
        pattern = r'(\d+),([\d]+\.[\d]*),(\d+),([UMH]+)\s'
        hash_of_sort_detail = {}
        for i in range(start_round,end_round + 1):
            if i in excepts_rounds:
                continue
            hash_of_i = {}
            infile_path = input_file_dir + os.sep + "detail_" + str(i) + ".csv"
            infile = open(infile_path, "r")
            line = infile.readline()
            while line:
                match = re.search(pattern, line)
                thread_no = int(match.group(1))
                generation_step = float(match.group(2))
                seq = match.group(4)

                if generation_step not in hash_of_i.keys():
                    hash_of_i[generation_step] = {}
                hash_of_i[generation_step][thread_no] = seq

                line = infile.readline()
                if not line:
                    break
            hash_of_i_new = {}
            for (key, val) in hash_of_i.items():
                sorted_hash_of_gen = sorted(val.iteritems(), key=lambda d: d[0])
                hash_of_i_new[key] = sorted_hash_of_gen
            hash_of_sort_detail[i] = hash_of_i_new
            infile.close()
        generation_steps = hash_of_sort_detail[start_round].keys()

        for gen_step in generation_steps:
            gen_step_str=str(gen_step).replace(".","_")
            out_file_path = output_file_dir + os.sep + "gen_" + gen_step_str+ "_detail.csv"
            out_file = open(out_file_path, 'w')
            for round_i in range(start_round,end_round + 1):
                if round_i in excepts_rounds:
                    continue
                seq_list = hash_of_sort_detail[round_i][gen_step]
                for (index, item) in seq_list:
                    if index == 0:
                        out_file.write(str(round_i) + ",")
                    out_file.write(item)
                out_file.write("\n")
            out_file.close()
    def get_sites_index_arr_from_file(self,out_cpg_sites_origin_path): # 从bed文件中提取所有的位点位置形成一个数组
        cpg_indexs = []
        cpg_sites_file = open(out_cpg_sites_origin_path, "r")
        line = cpg_sites_file.readline()

        while line:
            index_line_arr = line.split(",")
            # 提取postion,写入position value \n.
            pos_index = int(index_line_arr[1])
            cpg_indexs.append(pos_index)
            line = cpg_sites_file.readline()
            if not line and line.strip() == "":
                break
        cpg_sites_file.close()

        return cpg_indexs

    def convert_sorted_result_to_bed(self,pos_list, detail_dir, bed_target_dir,gens=[]):# 从CpG_sites位点和值得原始文件结合输出的模拟的各位点状态统计最后模拟的各个generation的模拟甲基化水平,输出到相应的文件夹
        for gen_i in gens:
            print "now convert the  %f gen sorted result to bed!" % gen_i
            gen_step_str = str(gen_i).replace(".", "_")
            detail_file = open(detail_dir + os.sep + "gen_" + gen_step_str + "_detail.csv", "r")
            bed_file = open(bed_target_dir + os.sep + "gen_" +gen_step_str+ ".bed", "w")
            detail_line = detail_file.readline()
            index = 0
            methy_status = []
            while detail_line:
                index = index + 1
                line_arr = detail_line.split(",")
                seq = line_arr[1]

                for no, ch in enumerate(seq):
                    if ch == "M":
                        if index == 1:
                            methy_status.append(2)
                        else:
                            methy_status[no] = methy_status[no] + 2
                    elif ch == "H":
                        if index == 1:
                            methy_status.append(1)
                        else:
                            methy_status[no] = methy_status[no] + 1
                    elif ch == "U":
                        if index == 1:
                            methy_status.append(0)
                detail_line = detail_file.readline()
                if not detail_line:
                    break
            print "index is %d" % index
            methy_status_arr = np.array(methy_status)
            methy_mean_arr = methy_status_arr / (float(index) * 2)
            methy_mean_list = methy_mean_arr.tolist()

            len_of_mean_list = len(methy_mean_list)
            for mean_index in range(len_of_mean_list):
                line_to_wrt = str(pos_list[mean_index]) + " " + str(methy_mean_list[mean_index]) + "\n"
                bed_file.write(line_to_wrt)
            bed_file.close()
            detail_file.close()
            print "now finished convert the  %f gen sorted result to bed!" % gen_i
    def sort_to_bed(self,pos_list, detail_dir, bed_dir,gens):
        if not os.path.exists(bed_dir):
            os.makedirs(bed_dir)
        self.convert_sorted_result_to_bed(pos_list, detail_dir, bed_dir, gens)
    def calc_corr(self,bed_dir,rd_with_dir,rd_without_dir,gen_collection,d_max,calc_interval=False,ignore_d=False):
        if calc_interval==True:
            if not os.path.exists(rd_with_dir):
                os.makedirs(rd_with_dir)
        else:
            if not os.path.exists(rd_without_dir):
                os.makedirs(rd_without_dir)
        for item in gen_collection:
            print "now calc %s correlation!" % item
            item_str=str(item).replace(".","_")
            #输出文件
            out_R_d_without_interval_file_name=rd_without_dir+os.sep+"chr1_r_d_without_"+item_str+".csv"
            #输入的bed文件
            bed_input=bed_dir+os.sep+"gen_"+str(item)+".bed"
            #计算相关性,将d和rd输出到文件中,d:2-corr_end
            if calc_interval==True:
                out_R_d_with_interval_file_name=rd_with_dir+os.sep+"chr1_r_d_with_"+item_str+".csv"
                calc_correlation(bed_input,out_R_d_with_interval_file_name,d_max,True,ignore_d)
            if ignore_d==True:
                calc_correlation(bed_input,out_R_d_without_interval_file_name,d_max,True,True)
            else:
                calc_correlation(bed_input,out_R_d_without_interval_file_name,d_max,False,False)
            print "now finished calc %s correlation!" % item

    def read_bed_file_and_store_pos_to_a_struct(self,bedfile_path, ignore_d=False):
        struct_to_store = {}
        bed_file = open(bedfile_path, 'r')
        re_pattern = r'(\d+)\s([\d]+\.[\d]*)\s'
        line = bed_file.readline()
        index = 0
        while line:
            match = re.search(re_pattern, line)
            if match:
                index = index + 1
                if ignore_d == True:
                    pos = index
                else:
                    pos = int(match.group(1))
                methy_level = float(match.group(2))
                struct_to_store[pos] = methy_level
            line = bed_file.readline()
        bed_file.close()
        # struct_to_store=sorted(struct_to_store.items(), key=lambda d:d[0])
        return struct_to_store

    # 根据距离d筛选在hash表中存储的所有满足要求距离的CpG对,输出为[[a_1 a_2] [a_1' a_2']......]形式
    def filter_d_length_to_generate_CpG_pairs(self,CpG_pos_and_methy_struct, d):
        array_to_store_pairs = []
        for key in CpG_pos_and_methy_struct.keys():
            pos = key
            methy_level_1 = CpG_pos_and_methy_struct[key]
            pos_off_by_d = pos + d
            if CpG_pos_and_methy_struct.has_key(pos_off_by_d):
                methy_level_2 = CpG_pos_and_methy_struct[pos_off_by_d]
                array_to_store_pairs.append([methy_level_1, methy_level_2])
        return array_to_store_pairs

    # 20160403修改版,d只统计两个cpg位点间没有其他位点,并且两位点距离为d的cpg_pairs
    def filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(self,od, d, length_of_od, od_keys, od_vals):
        # 排好序的dict
        array_to_store_pairs = []
        # 只遍历到倒数第二个元素,否则+1就会出界
        for i in range(0, length_of_od - 1):
            # 后一个位点与前一个位点之间的距离为d
            pre_od_key = od_keys[i]
            post_od_key = od_keys[i + 1]
            # print "now d:%d,now i%d of len:%d" %(d,i,length_of_od)
            if pre_od_key + d == post_od_key:
                methy_level_1 = od_vals[i]
                methy_level_2 = od_vals[i + 1]
                array_to_store_pairs.append([methy_level_1, methy_level_2])
        return array_to_store_pairs

    # 根据第二种方法计算相关系数r(d)
    def calc_C_d_by_pearson_correlation(self,CpG_pairs):
        sum_pre = 0.0
        sum_post = 0.0
        for pair in CpG_pairs:
            sum_pre = sum_pre + pair[0]
            sum_post = sum_post + pair[1]
        length = len(CpG_pairs)
        mean_1 = sum_pre / length
        mean_2 = sum_post / length

        sum_up = 0.0
        sum_down_left = 0.0
        sum_down_right = 0.0
        for pair in CpG_pairs:
            X_i = pair[0]
            Y_i = pair[1]
            sum_up = sum_up + (X_i - mean_1) * (Y_i - mean_2)
            sum_down_left = sum_down_left + (X_i - mean_1) * (X_i - mean_1)
            sum_down_right = sum_down_right + (Y_i - mean_2) * (Y_i - mean_2)
        sum_down = math.sqrt(sum_down_left * sum_down_right)
        if sum_down == 0:
            return -1
        r_d = sum_up / sum_down
        return r_d
    def calc_correlation(self,bed_file_path, out_R_d_file_path, d_max, is_inter_with_other_cpg, ignore_d=False):
        CpG_pos_and_methy_struct = self.read_bed_file_and_store_pos_to_a_struct(bed_file_path, ignore_d)
        # 要计算的d的范围
        d_list = range(2, d_max)
        out_R_d_file = open(out_R_d_file_path, 'w')
        sorted_struct = collections.OrderedDict(sorted(CpG_pos_and_methy_struct.items()))
        od_keys = sorted_struct.keys()
        od_vals = sorted_struct.values()
        length_of_od = len(sorted_struct)

        for d in d_list:
            if is_inter_with_other_cpg:
                CpG_pairs = self.filter_d_length_to_generate_CpG_pairs(CpG_pos_and_methy_struct, d)
            else:
                CpG_pairs = self.filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(sorted_struct, d, length_of_od,
                                                                                           od_keys, od_vals)
            # out_file("cpg_pairs.csv",CpG_pairs)
            d_count = len(CpG_pairs)
            print "finish chr%s d=%d run , d_count=%d" % ("1", d, d_count)
            if len(CpG_pairs) == 0:
                print "passed d=%d" % d
                continue
            r_d = self.calc_C_d_by_pearson_correlation(CpG_pairs)
            if r_d != -1:
                line2 = str(d) + "," + str(r_d) + "\n"
                out_R_d_file.write(line2)
                print "finish chr%s d=%d run , r_d=%f" % ("1", d, r_d)
        out_R_d_file.close()

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
    #simulator.run()

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
    list_gen = [simulator.rounds[len(simulator.rounds) - 1] + 0.01]
    str_list_gen = []
    for item in list_gen:
        str_list_gen.append(str(item))

    simulator.calc_corr(bed_files_dir,rd_with_dir,rd_without_dir,str_list_gen,calc_d_max,calc_interval=calc_interval,ignore_d=ignore_d)
def random_col_simulation(function_util):
    random_propensity_list = function_util.set_collaborative_params(U_plus_in=0.0005,H_plus_in=0.0005,M_minus_in=0.07,H_minus_in=0.069 ,H_p_H_in=0.22,H_p_M_in=0.22,U_p_M_in=0.23 ,H_m_U_in=0.072,M_m_U_in=0.06)

    sim_rounds = range(1, 2)

    random_index = "random_simulation"
    RANDOM_OUTPUT_DIR = "data" + os.sep + random_index

    max_cpg_sites = 1000
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

    init_cell = function_util.generate_CpG_in_methylation_percent_UHM(max_cpg_sites, m_ratio,u_ratio)  # generate a cpg chain which have the methylation status
    simulator = Simulator(random_propensity_list, rounds=sim_rounds, out_dir=RANDOM_OUTPUT_DIR,
                      max_cpg_sites=max_cpg_sites, generations=generations, pos_list=pos_list,
                      multi_threads=multi_threads, init_cell=init_cell, nearby=nearby, max_cells=max_cells,
                      index=random_index, detail_for_timestep=detail_for_timestep)
    simulator.run()

    sorted_ratio_dir_name = "sorted_ratio"
    sorted_ratio_bk_dir_name = "sorted_ratio_bk"
    sorted_detail_dir_name = "sorted_detail"

    sorted_ratio_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_ratio_dir_name
    sorted_ratio_bk_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_ratio_bk_dir_name
    sort_detail_dir = RANDOM_OUTPUT_DIR + os.sep + sorted_detail_dir_name

    simulator.sort_the_simulaiton_result(RANDOM_OUTPUT_DIR, sorted_ratio_dir, sort_detail_dir, simulator.rounds[0],
                                         simulator.rounds[len(simulator.rounds) - 1], [], False)

    list_gen=[simulator.rounds[len(simulator.rounds) - 1] + 0.01]
    str_list_gen=[]
    for item in list_gen:
        str_list_gen.append(str(list_gen))
if __name__ == '__main__':
    function_util = FunctionUtil()
    traditional_simulation(function_util)
    #random_col_simulation(function_util)