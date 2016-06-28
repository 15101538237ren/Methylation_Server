# -*- coding:utf-8 -*-
import random,math
import numpy as np
import matplotlib.pyplot as plt
class FuncionUtil(object):
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
if __name__ == '__main__':
    #create function utility object
    function_util=FuncionUtil()

    #set param for create a pos list obey geometric distribution
    n_cpg_sites=500
    geometric_p=0.3
    plot=True
    cpg_max_pos, cpg_pos_list=function_util.construct_n_cpg_sites_for_exp_distribution(n_cpg_sites,geometric_p,plot=plot)