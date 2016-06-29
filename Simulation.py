# -*- coding:utf-8 -*-
import random,math,os,threading,datetime
import numpy as np

BED_FILE_PATH="chr1.bed"
THRESHOLD_U=0.2
THRESHOLD_M=0.8
N_STEP=100

REACTION_STATUS_HASH={0:'H',1:'M',2:'H',3:'U',4:'M',5:'M',6:'H',7:'U',8:'H'}
RIGHT_STATUS_OF_REACTION={0:"U",1:"H",2:"M",3:"H",4:"H",5:"H",6:"U",7:"H",8:"M"}
STATE_OF_COLLABOR_REACTION={4:"H",5:"M",6:"M",7:"U",8:"U"}
BASE_RATE_HASH={4:1,5:1,6:0,7:3,8:2} #collaboration reaction base reaction rate index_hash
COLLABORATION_INDEX=4
class Simulatior(object):
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
            self.num_sites_per_thread=1000
            self.max_threads = 250
            self.num_of_threads = int(math.ceil(float(len(init_cell)) / self.num_sites_per_thread)) # thread counts
            self.num_turns = int(math.ceil(float(self.num_of_threads) / self.max_threads)) # when the max_thread is not enough for one round simulation, split it in to num_turns
            self.threads = [] # thread array
            print "%d thread is needed!" % self.num_of_threads
    def run(self):
        starttime_new=datetime.datetime.now()
        for i_round in self.rounds:
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
                        t.start()
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
                        CpG_str_list=list(cell_collection[idx])
                        CpG_str_list[target_reaction_CpG_site]=REACTION_STATUS_HASH[reaction_id]
                        cell_collection[idx]=''.join(CpG_str_list)
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