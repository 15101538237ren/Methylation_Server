# -*- coding:utf-8 -*-
import random,re,math
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
class FunctionUtil(object):
    '''
    This is a function utilities class
    '''

    #Construction function
    def __init__(self):
        pass

    #plot d distribution in histogram: d is a list, bins is the number of the bin in the hist.
    def plot_d_distribution(self,d,bins=100):
        plt.hist(d, bins=bins)  # plt.hist passes it's arguments to np.histogram
        plt.title("Histogram with "+str(bins)+" bins")
        plt.show()

    #Construct a CpG chain which the distance distribution is according
    #to the exponent distribution, n: site count, p: exponet distribution param, return the max cpg position and the cpg position list
    def construct_n_cpg_sites_for_exp_distribution(self,n,p,plot=False):
        cpg_pos_list=[]        # store the cpg position list which follow the exponent distribution
        cpg_pos=random.randint(1,10)  #give a init position of the first cpg site
        cpg_pos_list.append(cpg_pos)

        d=(np.random.geometric(p=p, size=n)+1).tolist()
        if plot:
            self.plot_d_distribution(d)
        for i in range(n-1):
            cpg_pos = d[i] + cpg_pos
            cpg_pos_list.append(cpg_pos)
        return cpg_pos,cpg_pos_list
    def get_nsite_for_d_distance(self,n,d):
        cpg_pos_list=[]        # store the cpg position list which follow the exponent distribution
        pos_now=1
        for i in range(n):
            cpg_pos_list.append(pos_now)
            pos_now=pos_now+d
        pos_now=pos_now-d
        return pos_now,cpg_pos_list
    # get all positions forming a list from the bed file given bed file path and the max CpG sites to read
    def get_pos_list_from_bed_file(self,bed_file_path,max_cpg_sites=0):
        methy_file = open(bed_file_path,'r')
        pos_list=[]
        line = methy_file.readline()
        methy_patern=r'(\d+)(\s)([\d]+\.[\d]*)'
        index=0
        print "now reading bed file!"
        while line:
            if max_cpg_sites>0:
                if index >= max_cpg_sites:
                    break
            match = re.search(methy_patern,line)
            methy_pos=int(match.group(1))
            line = methy_file.readline()
            pos_list.append(methy_pos)
            index=index+1
            if not line:
                break
        methy_file.close()
        return index,pos_list
    # generate the cpg site string which follow the m_ratio and u_ratio
    def generate_CpG_in_methylation_percent_UHM(self,CpG_sites_counts, m_ratio, u_ratio):
        CpG_str = "" # null string

        for i in range(1, CpG_sites_counts + 1): # random a number which determin the status of the target site
            random_num = random.random()
            if random_num <= m_ratio:
                CpG_str += "M"
            elif random_num > m_ratio and random_num <= m_ratio + u_ratio:
                CpG_str += "U"
            else:
                CpG_str += "H"
        # 返回生成的状态字符串,长度=CpG_sites_counts
        return CpG_str
    def set_standard_params_from(self,**reaction_hash):
        U_plus_in=float(reaction_hash.get("U_plus_in",0.0))
        H_plus_in=float(reaction_hash.get("H_plus_in",0.0))
        M_minus_in=float(reaction_hash.get("M_minus_in",0.0))
        H_minus_in=float(reaction_hash.get("H_minus_in",0.0))
        return [U_plus_in,H_plus_in,M_minus_in,H_minus_in]
    #set the standard reaction's parameters and return the parameter list
    def set_standard_params(self,U_plus_in,H_plus_in,M_minus_in,H_minus_in):
        return [U_plus_in,H_plus_in,M_minus_in,H_minus_in]
    #set the collaborative reaction's parameters given the reaction param hash and return the parameter list
    def set_collaborative_params(self,**reaction_hash):
        U_plus_in=float(reaction_hash.get("U_plus_in",0.0))
        H_plus_in=float(reaction_hash.get("H_plus_in",0.0))
        M_minus_in=float(reaction_hash.get("M_minus_in",0.0))
        H_minus_in=float(reaction_hash.get("H_minus_in",0.0))
        H_p_H_in=float(reaction_hash.get("H_p_H_in",0.0))
        H_p_M_in=float(reaction_hash.get("H_p_M_in",0.0))
        U_p_M_in=float(reaction_hash.get("U_p_M_in",0.0))
        H_m_U_in=float(reaction_hash.get("H_m_U_in",0.0))
        M_m_U_in=float(reaction_hash.get("M_m_U_in",0.0))
        return [U_plus_in, H_plus_in, M_minus_in, H_minus_in, H_p_H_in, H_p_M_in, U_p_M_in, H_m_U_in, M_m_U_in]
    def get_u_h_m_portion(self,bed_file_path,max_cpg_sites=0,u_thresh=0.2,m_thresh=0.8):
        methy_file = open(bed_file_path,'r')
        line = methy_file.readline()
        methy_patern=r'(\d+)(\s)([\d]+\.[\d]*)'
        index=0
        print "now reading bed file!"
        u_count=0
        m_count=0
        h_count=0
        while line:
            if max_cpg_sites>0:
                if index >= max_cpg_sites:
                    break
            if index % 10000 == 0:
                print "now readed %d lines!" % index
            match = re.search(methy_patern,line)
            methy_level=float(match.group(3))
            if methy_level<u_thresh:
                u_count=u_count+1
            elif methy_level>m_thresh:
                m_count=m_count+1
            else:
                h_count=h_count+1
            line = methy_file.readline()
            index=index+1
            if not line:
                break
        methy_file.close()
        u_ratio=float(u_count)/float(index)
        h_ratio=float(h_count)/float(index)
        m_ratio=float(m_count)/float(index)
        print "u: %f, h: %f, m: %f" %(u_ratio,h_ratio,m_ratio)
        return [u_ratio,h_ratio,m_ratio]
    def get_CGI_from(self,cgi_file_path,chr_no,start_pos=-1,end_pos=-1):
        #给定染色体编号,起始结束位点.得出某条染色体的所有CGI的数据,返回一个在start_pos和end_pos范围内的CG的List[[start1,end1],[start2,end2]...]
        cgi_file=open(cgi_file_path,"r")
        cgi_list=[]
        line=cgi_file.readline()
        methy_patern=r'(\d+)\schr([\d|X|Y]+)\s(\d+)\s(\d+)\sCpG:\s(\d+)'
        start_default=(start_pos < 0 or end_pos <=0 )

        while line:
            match = re.search(methy_patern,line)
            if match:
                chr_now=match.group(2)
                start=int(match.group(3))
                end=int(match.group(4))
                if chr_now==chr_no and ( start_default or (start>=start_pos and end <= end_pos)):
                    cpg_count=int(match.group(5))
                    cgi_list.append([start,end,cpg_count])
            line=cgi_file.readline()
        cgi_file.close()
        return cgi_list
    def sort_experiment_and_simulation_bed_by_cpg_methylation(self,cgi_list,exp_file_path,sim_file_path,cpg_counter_thresh,out_info_file_path,out_sorted_methylation_dir):
        if not os.path.exists(out_sorted_methylation_dir):
            os.makedirs(out_sorted_methylation_dir)
        #给定染色体编号,起始结束位点.得出某条染色体的所有CGI的数据,返回一个在start_pos和end_pos范围内的CG的List[[start1,end1],[start2,end2]...]
        exp_methylation_hash={}
        sim_methylation_hash={}
        methy_patern=r'(\d+)\s([\d]+\.[\d]*)'
        exp_bed_file=open(exp_file_path,"r")
        exp_line = exp_bed_file.readline()

        sim_bed_file=open(sim_file_path,"r")
        sim_line = sim_bed_file.readline()
        index=0
        cgi_index=0
        cgi_list_end=False
        print "now reading bed file!"
        while exp_line and sim_line:
            if index % 10000 == 0:
                print "readed %d line!" % index
            exp_match = re.search(methy_patern,exp_line)
            sim_match = re.search(methy_patern,sim_line)
            if exp_match and sim_match:
                exp_methy_pos=int(exp_match.group(1))
                sim_methy_pos=int(sim_match.group(1))
                if exp_methy_pos!=sim_methy_pos:
                    print "exp_methy_pos is not equal to sim_methy_pos!"
                    break
                else:
                    exp_methy_level=float(exp_match.group(2))
                    sim_methy_level=float(sim_match.group(2))
                    while exp_methy_pos > cgi_list[cgi_index][1]:
                        # now methylation position >CGI right bound
                        cgi_index=cgi_index+1
                        if cgi_index>=len(cgi_list):
                            cgi_list_end=True
                            break
                    if cgi_list_end:
                        break
                    if exp_methy_pos >= cgi_list[cgi_index][0]:
                        if cgi_index not in exp_methylation_hash.keys():
                            exp_methylation_hash[cgi_index]=[]
                            sim_methylation_hash[cgi_index]=[]
                        exp_methylation_hash[cgi_index].append([exp_methy_pos,exp_methy_level])
                        sim_methylation_hash[cgi_index].append([exp_methy_pos,sim_methy_level])
            exp_line = exp_bed_file.readline()
            sim_line = sim_bed_file.readline()
            index=index+1
        exp_bed_file.close()
        sim_bed_file.close()

        # for key,val in exp_methylation_hash.items():
        #     list_len=len(val)
        #     if list_len <= cpg_counter_thresh:
        #         del exp_methylation_hash[key]
        #         del sim_methylation_hash[key]

        sort_sim_hash={}
        sort_exp_hash={}
        for key,val in sim_methylation_hash.items():
            list_len=len(val)
            sum_of_sim_list=0.0
            sum_of_exp_list=0.0
            for i in range(list_len):
                sum_of_sim_list=sum_of_sim_list+val[i][1]
                sum_of_exp_list=sum_of_exp_list+exp_methylation_hash[key][i][1]
            mean_of_sim_list=sum_of_sim_list/float(list_len)
            mean_of_exp_list=sum_of_exp_list/float(list_len)
            sort_sim_hash[key]=mean_of_sim_list
            sort_exp_hash[key]=mean_of_exp_list
        sorted_sim_dict= sorted(sort_sim_hash.iteritems(), key=lambda d:d[1], reverse = True)
        out_file=open(out_info_file_path,"w")
        head_line="CGI_id,CGI start site,CGI end site,mean experiment methylation level,mean simulation methylation level,CpG count of the CGI, experiment detected CpG count\n"
        out_file.write(head_line)

        for i in range(len(sorted_sim_dict)):
            key,val=sorted_sim_dict[i]

            cgi_start=cgi_list[key][0]
            cgi_end=cgi_list[key][1]
            total_cgi=cgi_list[key][2]
            detected_cgi=len(exp_methylation_hash[key])
            line_to_wrt=str(key)+","+str(cgi_start)+","+str(cgi_end)+","+str(sort_exp_hash[key])+","+str(val)+","+str(total_cgi)+","+str(detected_cgi)+"\n"
            out_file.write(line_to_wrt)
            out_exp_sim_file_path=out_sorted_methylation_dir+os.sep+str(i)+"_"+str(key)+".csv"
            out_exp_sim_file=open(out_exp_sim_file_path,"w")
            # head_of_rank_i="methylation position,experiment methylation level,simulation methylation level\n"
            head_of_rank_i="column,row,methy,label\n"
            out_exp_sim_file.write(head_of_rank_i)
            rows=4
            cols=int(math.ceil(float(len(exp_methylation_hash[key]))/float(rows)))
            for j in range(len(exp_methylation_hash[key])):
                #methy_pos=exp_methylation_hash[key][j][0]
                row_no=rows-(int(j/cols))-1
                col_no=j % cols
                exp_methy_lev=exp_methylation_hash[key][j][1]
                if exp_methy_lev > 0.8:
                    status="Methylated"
                elif exp_methy_lev > 0.2:
                    status="Half-methylated"
                else:
                    status="Unmethylated"
                line_of_rank_i_to_wrt=str(col_no)+","+str(row_no)+","+status+","+"Experiment"+"\n"
                out_exp_sim_file.write(line_of_rank_i_to_wrt)
            for k in range(len(sim_methylation_hash[key])):
                sim_methy_lev=sim_methylation_hash[key][k][1]
                #line_of_rank_i_to_wrt=str(methy_pos)+","+str(exp_methy_lev)+","+str(sim_methy_lev)+"\n"
                row_no=rows-(int(k/cols))-1
                col_no=k % cols
                if sim_methy_lev > 0.8:
                    status="Methylated"
                elif sim_methy_lev > 0.2:
                    status="Half-methylated"
                else:
                    status="Unmethylated"
                line_of_rank_i_to_wrt=str(col_no)+","+str(row_no)+","+status+","+"Simulation"+"\n"
                out_exp_sim_file.write(line_of_rank_i_to_wrt)

            out_exp_sim_file.close()
        out_file.close()
    def sort_CGI(self,bed_file_path,out_file_path,cgi_list,max_cpg_sites=0):
        methylation_hash={}
        bed_file=open(bed_file_path,"r")
        line = bed_file.readline()
        methy_patern=r'(\d+)\s([\d]+\.[\d]*)'
        index=0
        cgi_index=0
        cgi_list_end=False
        print "now reading bed file!"
        while line:
            if max_cpg_sites>0:
                if index >= max_cpg_sites:
                    break
            match = re.search(methy_patern,line)
            if match:
                methy_pos=int(match.group(1))
                methy_level=float(match.group(2))
                try:
                    while methy_pos > cgi_list[cgi_index][1]:
                        # now methylation position >CGI right bound
                        cgi_index=cgi_index+1
                        if cgi_index>=len(cgi_list):
                            cgi_list_end=True
                            break
                except IndexError,e:
                    print e
                if cgi_list_end:
                    break
                if methy_pos >= cgi_list[cgi_index][0]:
                    if cgi_index not in methylation_hash.keys():
                        methylation_hash[cgi_index]=[]
                    methylation_hash[cgi_index].append([methy_pos,methy_level])
            line = bed_file.readline()

            index=index+1
        bed_file.close()
        #sort the result CGI by their mean methylation level
        sort_hash={}
        for key,val in methylation_hash.items():
            list_len=len(val)
            sum_of_list=0.0
            for i in range(list_len):
                sum_of_list=sum_of_list+val[i][1]
            mean_of_list=sum_of_list/float(list_len)
            sort_hash[key]=mean_of_list
        sorted_dict= sorted(sort_hash.iteritems(), key=lambda d:d[1], reverse = True)

        out_file=open(out_file_path,"w")
        # head_line="CGI start site,CGI end site,mean methylation level,CpG count of the CGI, experiment detected CpG count\n"
        # out_file.write(head_line)
        for i in range(len(sorted_dict)):
            key,val=sorted_dict[i]
            cgi_start=cgi_list[key][0]
            cgi_end=cgi_list[key][1]
            total_cgi=cgi_list[key][2]
            detected_cgi=len(methylation_hash[key])
            line_to_wrt=str(cgi_start)+","+str(cgi_end)+","+str(val)+","+str(total_cgi)+","+str(detected_cgi)+"\n"
            out_file.write(line_to_wrt)
        out_file.close()
        print "now finished sort CGI!"
        return sorted_dict
    def sort_CGI_by_sorted_dict(self,bed_file_path,out_file_path,cgi_list,sorted_site_list,max_cpg_sites=0):
        methylation_hash={}
        bed_file=open(bed_file_path,"r")
        line = bed_file.readline()
        methy_patern=r'(\d+)\s([\d]+\.[\d]*)'
        index=0
        cgi_index=0
        cgi_list_end=False
        print "now reading bed file!"
        while line:
            if max_cpg_sites>0:
                if index >= max_cpg_sites:
                    break
            match = re.search(methy_patern,line)
            if match:
                methy_pos=int(match.group(1))
                methy_level=float(match.group(2))
                try:
                    while methy_pos > cgi_list[cgi_index][1]:
                        # now methylation position >CGI right bound
                        cgi_index=cgi_index+1
                        if cgi_index>=len(cgi_list):
                            cgi_list_end=True
                            break
                except IndexError,e:
                    print e
                if cgi_list_end:
                    break
                if methy_pos >= cgi_list[cgi_index][0]:
                    if cgi_index not in methylation_hash.keys():
                        methylation_hash[cgi_index]=[]
                    methylation_hash[cgi_index].append([methy_pos,methy_level])
            line = bed_file.readline()

            index=index+1
        bed_file.close()
        #sort the result CGI by their mean methylation level
        sort_hash={}
        for key,val in methylation_hash.items():
            list_len=len(val)
            sum_of_list=0.0
            for i in range(list_len):
                sum_of_list=sum_of_list+val[i][1]
            mean_of_list=sum_of_list/float(list_len)
            sort_hash[key]=mean_of_list

        out_file=open(out_file_path,"w")
        # head_line="CGI start site,CGI end site,mean methylation level,CpG count of the CGI, experiment detected CpG count\n"
        # out_file.write(head_line)
        methylation_level_list=[]
        for i in range(len(sorted_site_list)):
            key=sorted_site_list[i]
            cgi_start=cgi_list[key][0]
            cgi_end=cgi_list[key][1]
            total_cgi=cgi_list[key][2]
            if key in sort_hash.keys():
                detected_cgi=len(methylation_hash[key])
                line_to_wrt=str(cgi_start)+","+str(cgi_end)+","+str(sort_hash[key])+","+str(total_cgi)+","+str(detected_cgi)+"\n"
                out_file.write(line_to_wrt)
                methylation_level_list.append(sort_hash[key])
            else:
                #没有跟2细胞期的结构匹配上的CGI,如其他细胞期检测到,但是2细胞期没有检测到的
                line_to_wrt=str(cgi_start)+","+str(cgi_end)+","+str(-1.0)+","+str(total_cgi)+","+str(-1)+"\n"
                out_file.write(line_to_wrt)
                methylation_level_list.append("NA")
        out_file.close()
        print "now finished sort CGI!"
        return methylation_level_list
def sort_CGI_by_2cell():
    rootdir = "/Users/Ren/PycharmProjects/Methylation_Server/"
    cgi_file_path="cpgIslandExt.txt"
    out_dir="CGI_sort_data_without_head"
    if not os.path.exists(rootdir+os.sep+out_dir):
        os.makedirs(rootdir+os.sep+out_dir)
    function_util=FunctionUtil()
    cgi_list=function_util.get_CGI_from(cgi_file_path,"1")

    cell_2_path="2 cell paternal"
    bed_file_path=cell_2_path+os.sep+"sp_bed"+os.sep+"chr1.bed"
    out_file_path=out_dir+os.sep+cell_2_path+".csv"
    print cell_2_path
    sorted_dict=function_util.sort_CGI(bed_file_path,out_file_path,cgi_list)
    sorted_site_list=[]
    methylation_level={}
    methylation_level_tmp_list=[]
    for i in range(len(sorted_dict)):
        key,val=sorted_dict[i]
        sorted_site_list.append(key)
        methylation_level_tmp_list.append(val)
    methylation_level[cell_2_path]=methylation_level_tmp_list
    len_of_the_result=len(methylation_level_tmp_list)
    input_dir=rootdir+os.sep+"matlab/mouse_bed"
    for parent,dirnames,filenames in os.walk(input_dir):    #三个参数：分别返回1.父目录 2.所有文件夹名字（不含路径） 3.所有文件名字
        for dirname in  dirnames:                       #输出文件夹信息
            if dirname =="sp_bed":
                last_names=parent.split("/")
                parent_dir=last_names[len(last_names)-1]
                print parent_dir
                bed_file_path=parent+os.sep+dirname+os.sep+"chr1.bed"
                out_file_path=out_dir+os.sep+parent_dir+".csv"
                methylation_level_tmp_list=function_util.sort_CGI_by_sorted_dict(bed_file_path,out_file_path,cgi_list,sorted_site_list)
                methylation_level[parent_dir]=methylation_level_tmp_list
                print "finishe handling %s" % parent_dir
    # out_file_path="cells_methy_by_2cell_rank.csv"
    # out_file=open(out_file_path,"w")
    # dir_name_list=["2 cell paternal","2 cell maternal","4 cell paternal","4 cell maternal","E6.5 paternal","E6.5 maternal","E7.5 paternal","E7.5 maternal","ICM maternal","PCG E13.5 male","PCG E13.5 female","sperm","oocyte"]
    # for i in range(len_of_the_result):
    #     if i==0:
    #         str_to_wrt=",".join(dir_name_list)
    #         out_file.write(str_to_wrt+"\n")
    #     tmp_arr=[]
    #     for key in dir_name_list:
    #         tmp_arr.append(str(methylation_level[key][i]))
    #     str_to_wrt=",".join(tmp_arr)
    #     out_file.write(str_to_wrt+"\n")
    # out_file.close()

def sort_all_result():
    #create function utility object
    rootdir = "/Users/Ren/PycharmProjects/Methylation_Server/"
    cgi_file_path="cpgIslandExt.txt"
    out_dir="CGI_sort_data"
    if not os.path.exists(rootdir+os.sep+out_dir):
        os.makedirs(rootdir+os.sep+out_dir)
    input_dir=rootdir+os.sep+"matlab/mouse_bed"

    function_util=FunctionUtil()
    cgi_list=function_util.get_CGI_from(cgi_file_path,"1")
    for parent,dirnames,filenames in os.walk(input_dir):    #三个参数：分别返回1.父目录 2.所有文件夹名字（不含路径） 3.所有文件名字
        for dirname in  dirnames:                       #输出文件夹信息
            if dirname =="sp_bed":
                last_names=parent.split("/")
                parent_dir=last_names[len(last_names)-1]
                print parent_dir
                bed_file_path=parent+os.sep+dirname+os.sep+"chr1.bed"
                out_file_path=out_dir+os.sep+parent_dir+"_CGI_result.csv"

                function_util.sort_CGI(bed_file_path,out_file_path,cgi_list)
                print "finishe handling %s" % parent_dir
# def read_exp_and_sim_generate_r_dataframe(csv_file_path,out_file_path):
#     sim_bed_file=open(sim_file_path,"r")
#     sim_line = sim_bed_file.readline()
if __name__ == '__main__':
    # sort_all_result()
    #sort_CGI_by_2cell()
    cgi_file_path="cpgIslandExt.txt"
    function_util=FunctionUtil()
    cgi_list=function_util.get_CGI_from(cgi_file_path,"1")
    indir="E75_paternal_compare"
    exp_bed_path="Mouse_experiment_bed"+os.sep+"E75_mc_CG_paternal_chr1.bed"
    sim_bed_path=indir+os.sep+"E75_paternal_gen_49_95.bed"
    out_info_path=indir+os.sep+"exp_and_sim_CGI_sort_result.csv"
    out_exp_sim_dir=indir+os.sep+"exp_and_sim_CGI"
    function_util.sort_experiment_and_simulation_bed_by_cpg_methylation(cgi_list,exp_bed_path,sim_bed_path,10,out_info_path,out_exp_sim_dir)

    # function_util.get_u_h_m_portion(bed_file_path,u_thresh=0.2,m_thresh=0.8)
    #
    # #set param for create a pos list obey geometric distribution
    # n_cpg_sites=500
    # geometric_p=0.3
    # plot=True
    # cpg_max_pos, cpg_pos_list=function_util.construct_n_cpg_sites_for_exp_distribution(n_cpg_sites,geometric_p,plot=plot)
    #
    # #the site origin ratio
    # m_ratio=0.181214
    # h_ratio=0.427782
    # u_ratio=0.391004
    #
    # #generate a cpg chain which have the methylation status
    # CpG_str=function_util.generate_CpG_in_methylation_percent_UHM(n_cpg_sites,m_ratio,u_ratio)
    # print CpG_str

