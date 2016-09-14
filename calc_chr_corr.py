# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle
import numpy as np
def read_bed_file_and_store_pos_to_a_struct(bedfile_path, ignore_d=False):
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
def filter_d_length_to_generate_CpG_pairs(CpG_pos_and_methy_struct, d): # 根据距离d筛选在hash表中存储的所有满足要求距离的CpG对,输出为[[a_1 a_2] [a_1' a_2']......]形式
        array_to_store_pairs = []
        for key in CpG_pos_and_methy_struct.keys():
            pos = key
            methy_level_1 = CpG_pos_and_methy_struct[key]
            pos_off_by_d = pos + d
            if CpG_pos_and_methy_struct.has_key(pos_off_by_d):
                methy_level_2 = CpG_pos_and_methy_struct[pos_off_by_d]
                array_to_store_pairs.append([methy_level_1, methy_level_2])
        return array_to_store_pairs
def filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(d, length_of_od, od_keys, od_vals): # 20160403修改版,d只统计两个cpg位点间没有其他位点,并且两位点距离为d的cpg_pairs
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
def calc_C_d_by_pearson_correlation(CpG_pairs): # 根据第二种方法计算相关系数r(d)
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
#根据第一种方法计算C(d)的值
def calc_C_d_by_prob_N1_N2(CpG_pairs,threshold_to_distinct_=0.5):
    N1=N2=N3=N4=0
    for pair in CpG_pairs:
        first=0
        second=0
        if pair[0]>=0.5:
            first=1
        if pair[1]>=0.5:
            second=1
        #N1情况
        if first==0 and second==0:
            N1=N1+1.0
        elif first==0 and second==1:
            N2=N2+1.0
        elif first==1 and second==0:
            N3=N3+1.0
        elif first==1 and second==1:
            N4=N4+1.0
    if (N1+N2)==0 or (N3+N4)==0:
        return -1
    C_d=(N1/(N1+N2)+N4/(N3+N4))-1.0

    return C_d
def calc_correlation(chr_no,bed_file_path, out_R_d_file_path, d_max, is_inter_with_other_cpg,rd=1, ignore_d=False):
        CpG_pos_and_methy_struct = read_bed_file_and_store_pos_to_a_struct(bed_file_path, ignore_d)
        # 要计算的d的范围
        d_list = range(1, d_max)
        out_R_d_file = open(out_R_d_file_path, 'w')
        sorted_struct = collections.OrderedDict(sorted(CpG_pos_and_methy_struct.items()))
        od_keys = sorted_struct.keys()
        od_vals = sorted_struct.values()
        length_of_od = len(sorted_struct)

        for d in d_list:
            if is_inter_with_other_cpg:
                CpG_pairs = filter_d_length_to_generate_CpG_pairs(CpG_pos_and_methy_struct, d)
            else:
                CpG_pairs = filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(d, length_of_od,
                                                                                           od_keys, od_vals)
            d_count = len(CpG_pairs)
            print "finish chr%s d=%d run , d_count=%d" % (str(chr_no), d, d_count)
            if len(CpG_pairs) == 0:
                print "passed d=%d" % d
                continue
            if rd:
                r_d = calc_C_d_by_pearson_correlation(CpG_pairs)
            else:
                r_d = calc_C_d_by_prob_N1_N2(CpG_pairs)
            if r_d != -1:
                line2 = str(d) + "," + str(r_d) + "\n"
                out_R_d_file.write(line2)
                print "finish chr%s d=%d run , r_d=%f" % (str(chr_no), d, r_d)
        out_R_d_file.close()
if __name__ == '__main__':
    input_dir = "human_splitted_bed"
    chr_no_list = range(1, 2)
    d_max = 1000
    is_inter_with_other_cpg = True
    ignore_d = True
    standard=False
    standard_file_end=["_plus","_minus"]
    out_corr_dir="mouse"+os.sep+"corr"
    # 如果文件夹不存在则应先创建
    if not os.path.exists(out_corr_dir):
        os.makedirs(out_corr_dir)
    for chr_i in chr_no_list:
        if standard:
            for file_end in standard_file_end:
                bed_file=input_dir+os.sep+"chr"+str(chr_i)+file_end+".bed"
                our_rd_path=out_corr_dir+os.sep+"chr"+str(chr_i)+file_end+"_without.csv"
                calc_correlation(chr_i,bed_file,our_rd_path,d_max,is_inter_with_other_cpg,ignore_d)
        else:
            bed_file = "chr1.bed"#input_dir + os.sep + "chr" + str(chr_i) + ".bed"
            our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without_cd_intuitive.csv"
            calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,rd=1, ignore_d=ignore_d)