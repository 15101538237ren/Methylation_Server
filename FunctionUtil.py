# -*- coding:utf-8 -*-
import random,re
import numpy as np
import matplotlib.pyplot as plt
class FunctionUtil(object):
    '''
    This is a function utilities class
    '''

    #Construction function
    def __init__(self):
        pass
    def plot_d_distribution(self,d,bins=100):
        plt.hist(d, bins=bins)  # plt.hist passes it's arguments to np.histogram
        plt.title("Histogram with "+str(bins)+" bins")
        plt.show()
    def construct_n_cpg_sites_for_exp_distribution(self,n,p,plot=False):
        #Construct a CpG chain which the distance distribution is according
        #to the exponent distribution, n: site count, lamda: exponet distribution param, return the max cpg position
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

    # generate the cpg site which follow the m_ratio and u_ratio
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
    def set_standard_params(self,U_plus_in,H_plus_in,M_minus_in,H_minus_in):
        return [U_plus_in,H_plus_in,M_minus_in,H_minus_in]
    def set_collaborative_params(self,**reaction_hash): # 设置协作反应的参数
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
    def get_pos_list_from_bed_file(self,bed_file_path,max_cpg_sites):
        methy_file = open(bed_file_path,'r')
        #格式[[position methy_val],[position methy_val]]
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

if __name__ == '__main__':
    #create function utility object
    function_util=FunctionUtil()

    #set param for create a pos list obey geometric distribution
    n_cpg_sites=500
    geometric_p=0.3
    plot=True
    cpg_max_pos, cpg_pos_list=function_util.construct_n_cpg_sites_for_exp_distribution(n_cpg_sites,geometric_p,plot=plot)

    #the site origin ratio
    m_ratio=0.181214
    h_ratio=0.427782
    u_ratio=0.391004

    #generate a cpg chain which have the methylation status
    CpG_str=function_util.generate_CpG_in_methylation_percent_UHM(n_cpg_sites,m_ratio,u_ratio)