# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle
import numpy as np

TYPE_CPG_TXT="CpG_txt"
TYPE_NORMAL="normal"
TYPE_BS_CG_TXT="TYPE_BS_CG_TXT"
TYPE_ENCODE_BED="TYPE_ENCODE_BED"
TYPE_WGBS="TYPE_WGBS"
TYPE_ALL_C="TYPE_ALL_C"
def variable_wig_data_extract_to_methy_normal(in_file_path,out_dir_path):
    if not os.path.exists(out_dir_path):
        os.makedirs(out_dir_path)
    #读入的文件路径
    raw_file = open(in_file_path,'r')

    #文件样例格式:chr1  3025349	3025350	0.6	3	2.染色体编号,CpG起始位点,CpG结束位点,甲基化水平,甲基化的reads数,未甲基化的reads数

    #写入的文件路径
    out_file_pre_path=out_dir_path+os.sep+"chr"

    methy_file1 = open(out_file_pre_path+"1.bed",'w')

    count=0
    chr_no=0
    pattern=r'(\d+)\s([\d]+\.[\d]*)'

    line=raw_file.readline()
    while line:
        if count % 1000000==0:
            print "%d lines of %s was processed" %(count,in_file_path)
        match_p = re.search(pattern,line)
        if not match_p:
            chr_pattern=r'variableStep chrom=chr(\d+)'
            match_chr=re.search(chr_pattern,line)
            if match_chr:
                print "now chr "+match_chr.group(1)
                out_file_path=out_file_pre_path+match_chr.group(1)+".bed"
                methy_file1.close()
                methy_file1=open(out_file_path,"w")
        else:
            count=count+1
            if count % 2 !=0:
                methy_pos = int(match_p.group(1))
                value = float(match_p.group(2))
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
        line=raw_file.readline()
    methy_file1.close()
    raw_file.close()
    print "finish %s chr%s data processing!" %(in_file_path,chr_no)

def bed_data_extract_to_methy(in_file_path,outdir_path,pattern_type):
    #如果文件夹不存在则应先创建
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)
    #读入的文件路径
    raw_file = open(in_file_path,'r')

    #文件样例格式:chr1  3025349	3025350	0.6	3	2.染色体编号,CpG起始位点,CpG结束位点,甲基化水平,甲基化的reads数,未甲基化的reads数

    #写入的文件路径
    #写入的文件路径
    out_file_pre_path=outdir_path+os.sep+"chr"

    methy_file1 = open(out_file_pre_path+"1.bed",'w')

    count=0

    chr_no=0
    chr_no_index=1
    type,pattern,methy_pos_index,value_index=pattern_generate(pattern_type)

    line=raw_file.readline()
    while line:
        count=count+1
        if count % 100000==0:
            print "%d lines of %s was processed" %(count,in_file_path)
        match_p = re.search(pattern,line)

        if match_p:
            now_chr_no=int(match_p.group(chr_no_index))
            if chr_no!=now_chr_no:
                chr_no=now_chr_no
                print "now chr "+str(now_chr_no)
                out_file_path=out_file_pre_path+str(now_chr_no)+".bed"
                methy_file1.close()
                methy_file1=open(out_file_path,"w")
            methy_pos = int(match_p.group(methy_pos_index))
            if type ==TYPE_CPG_TXT:
                if match_p.group(3)=="F":
                    value = float(match_p.group(value_index))/100.0
                    out_str = str(methy_pos) + "\t" + str(value) + "\n"
                    methy_file1.write(out_str)
            elif type==TYPE_BS_CG_TXT:
                methylated_read_count=int(match_p.group(value_index))
                total_read_count=int(match_p.group(value_index+1))
                if total_read_count!=0:
                    value = float(methylated_read_count)/float(total_read_count)
                else:
                    value = 0.0
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
            elif type==TYPE_ENCODE_BED or type==TYPE_WGBS:
                value = float(match_p.group(value_index))/100.0
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
            else:
                value = float(match_p.group(value_index))
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
        line=raw_file.readline()
    methy_file1.close()
    raw_file.close()
    print "finish %s chr%s data processing!" %(in_file_path,chr_no)
def pattern_generate(type):
    if type==TYPE_NORMAL:
        return normal_bed_pattern()
    elif type==TYPE_WGBS:
        return WGBS_bed_pattern()
    elif type==TYPE_ENCODE_BED:
        return ENCODE_bed_pattern()
    elif type==TYPE_BS_CG_TXT:
        return bs_cg_txt_pattern()
    elif type==TYPE_CPG_TXT:
        return CpG_txt_pattern()
def normal_bed_pattern():
    pattern=r'chr(\d+)\t(\d+)(\s)(\d+)(\s)([\d]+[\.[\d]*]?)(\s)(\d+)(\s)(\d+)'
    methy_pos_index=2
    value_index=6
    type=TYPE_NORMAL
    return type,pattern,methy_pos_index,value_index
def WGBS_bed_pattern():
    pattern=r'chr(\d+)\t(\d+)(\s)(\d+)(\s)([\d]+[\.[\d]*]?)(\s)(\d+)(\s)(\d+)'
    methy_pos_index=2
    value_index=6
    type=TYPE_WGBS
    return type,pattern,methy_pos_index,value_index
def ENCODE_bed_pattern():
    pattern=r'chr(\d+)\t(\d+)\t(\d+)\t\.\t(\d+)\t\+\t(\d+)\t(\d+)\t(\d+),(\d+),(\d+)\t(\d+)\t(\d+)'
    methy_pos_index=2
    value_index=11
    type=TYPE_ENCODE_BED
    return type,pattern,methy_pos_index,value_index
def bs_cg_txt_pattern():
    pattern=r'chr(\d+)\t(\d+)\s(\d+)\s(\d+)'
    methy_pos_index=2
    value_index=3
    type=TYPE_BS_CG_TXT
    return type,pattern,methy_pos_index,value_index
def CpG_txt_pattern():
    pattern=r'[^\s]+\s(\d+)\s(\d+)\s(F|R)\s(\d+)\s([\d]+[\.[\d]*]?)(\s)([\d]+[\.[\d]*]?)'
    methy_pos_index=2
    value_index=5
    type=TYPE_CPG_TXT
    return type,pattern,methy_pos_index,value_index
def bed_data_extract_to_methy_all_C(in_file_path,out_dir_path):
    if not os.path.exists(out_dir_path):
        os.makedirs(out_dir_path)
    #读入的文件路径
    raw_file = open(in_file_path,'r')

    #文件样例格式:chr1  3025349	3025350	0.6	3	2.染色体编号,CpG起始位点,CpG结束位点,甲基化水平,甲基化的reads数,未甲基化的reads数

    #写入的文件路径

    #写入的文件路径
    out_file_pre_path=out_dir_path+os.sep+"chr"

    methy_file1 = open(out_file_pre_path+"1.bed",'w')

    chr_no=0
    count=0

    pattern=r'(\d+)\t(\d+)\s([-|+])\s(CHG|CHH|CG)\s(\d+)\s(\d+)'

    line=raw_file.readline()
    methylation_total=0
    methylated_total=0
    methylation_position=0
    while line:
        count=count+1
        if count % 1000000==0:
            print "%d lines of %s was processed" %(count,in_file_path)
        match_p = re.search(pattern,line)

        if match_p:
            chr_no_of_line=int(match_p.group(1))
            if chr_no!=chr_no_of_line:
                chr_no=chr_no_of_line
                print "now chr "+str(chr_no_of_line)
                out_file_path=out_file_pre_path+str(chr_no_of_line)+".bed"
                methy_file1.close()
                methy_file1=open(out_file_path,"w")
            CG=match_p.group(4)
            if CG=="CG":
                strand=match_p.group(3)
                tmp_pos=int(match_p.group(2))
                tmp_methylated=int(match_p.group(5))
                tmp_total=int(match_p.group(6))
                #正链情况,下一次还得看看负链
                if strand=="+":
                    #上次只有+链CpG
                    if methylation_position!=0:
                        value=float(methylated_total)/float(methylation_total)
                        out_str = str(methylation_position) + "\t" + str(value) + "\n"
                        methy_file1.write(out_str)
                    methylation_position=tmp_pos
                    methylated_total=tmp_methylated
                    methylation_total=tmp_total
                elif strand=="-":
                    #上行有正链CpG
                    if methylation_position!=0:
                        methylated_total=methylated_total+tmp_methylated
                        methylation_total=methylation_total+tmp_total
                    #只有负链CpG
                    else:
                        methylation_position=tmp_pos
                        methylated_total=tmp_methylated
                        methylation_total=tmp_total
                    value=float(methylated_total)/float(methylation_total)
                    out_str = str(methylation_position) + "\t" + str(value) + "\n"
                    methy_file1.write(out_str)
                    methylation_total=0
                    methylated_total=0
                    methylation_position=0
        line=raw_file.readline()
    methy_file1.close()
    raw_file.close()
    print "finish %s chr%d data processing!" %(in_file_path,chr_no)

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
def calc_correlation(chr_no,bed_file_path, out_R_d_file_path, d_max, is_inter_with_other_cpg, ignore_d=False):
        CpG_pos_and_methy_struct = read_bed_file_and_store_pos_to_a_struct(bed_file_path, ignore_d)
        # 要计算的d的范围
        d_list = range(2, d_max)
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
            print "finish chr%d d=%d run , d_count=%d" % (chr_no, d, d_count)
            if len(CpG_pairs) == 0:
                print "passed d=%d" % d
                continue
            r_d = calc_C_d_by_pearson_correlation(CpG_pairs)
            if r_d != -1:
                line2 = str(d) + "," + str(r_d) + "\n"
                out_R_d_file.write(line2)
                print "finish chr%d d=%d run , r_d=%f" % (chr_no, d, r_d)
        out_R_d_file.close()
if __name__ == '__main__':
    d_max = 1000
    is_inter_with_other_cpg = False
    ignore_d = False


    # cpg_txt_base="Daniorerio"+os.sep+"GSE59916"
    # cpg_txt_list=os.listdir(cpg_txt_base)
    # #
    # for item in cpg_txt_list:
    #     if os.path.isfile(os.path.join(cpg_txt_base,item)):
    #         if item.split(".")[-1]=="txt":
    #             input_bed_path=cpg_txt_base+os.sep+item
    #             out_dir_path = item.split(".")[0]
    #             out_dir_full_path=cpg_txt_base+os.sep+out_dir_path
    #             print "now handle %s" % out_dir_path
    #             out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
    #             pattern_type=TYPE_CPG_TXT
    #             bed_data_extract_to_methy(input_bed_path,out_bed_path,pattern_type)
    #             out_corr_dir=out_dir_full_path+os.sep+"correlation"
    #             # 如果文件夹不存在则应先创建
    #             if not os.path.exists(out_corr_dir):
    #                 os.makedirs(out_corr_dir)
    #             chr_no_list=range(1,26)
    #             for chr_i in chr_no_list:
    #                 bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
    #                 our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
    #                 calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)


    # all_c_txt_base="Daniorerio"+os.sep+"GSE68087"
    # all_c_txt_list=os.listdir(all_c_txt_base)
    #
    # for item in all_c_txt_list:
    #     if os.path.isfile(os.path.join(all_c_txt_base,item)):
    #         if item.split(".")[-1]=="txt":
    #             input_bed_path=all_c_txt_base+os.sep+item
    #             chr_no_list = range(1,26)
    #             out_dir_path = item.split(".")[0]
    #             pattern_type=TYPE_ALL_C
    #             out_dir_full_path=all_c_txt_base+os.sep+out_dir_path
    #             if not os.path.exists(out_dir_full_path):
    #                 print "now handle %s" % out_dir_path
    #                 out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
    #                 bed_data_extract_to_methy_all_C(input_bed_path,out_bed_path)
    #                 out_corr_dir=out_dir_full_path+os.sep+"correlation"
    #                 # 如果文件夹不存在则应先创建
    #                 if not os.path.exists(out_corr_dir):
    #                     os.makedirs(out_corr_dir)
    #                 chr_no_list=range(1,26)
    #                 for chr_i in chr_no_list:
    #                     bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
    #                     our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
    #                     calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)

    # WGBS_type_base="Musmusculus"
    # WGBS_type_list=os.listdir(WGBS_type_base)
    #
    # for item in WGBS_type_list:
    #     if os.path.isfile(os.path.join(WGBS_type_base,item)):
    #         if item.split(".")[-1]=="txt":
    #             input_bed_path=WGBS_type_base+os.sep+item
    #             out_dir_path = item.split(".")[0]
    #             out_dir_full_path=WGBS_type_base+os.sep+out_dir_path
    #             print "now handle %s" % out_dir_path
    #             out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
    #             pattern_type=TYPE_WGBS
    #             bed_data_extract_to_methy(input_bed_path,out_bed_path,pattern_type)
    #             out_corr_dir=out_dir_full_path+os.sep+"correlation"
    #             # 如果文件夹不存在则应先创建
    #             if not os.path.exists(out_corr_dir):
    #                 os.makedirs(out_corr_dir)
    #             chr_no_list=range(1,20)
    #             for chr_i in chr_no_list:
    #                 bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
    #                 our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
    #                 calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)


    # bs_cg_txt_base="Homosapiens"+os.sep+"GSE70090"
    # bs_cg_txt_list=os.listdir(bs_cg_txt_base)
    # for item in bs_cg_txt_list:
    #     if os.path.isfile(os.path.join(bs_cg_txt_base,item)):
    #         if item.split(".")[-1]=="txt":
    #             input_bed_path=bs_cg_txt_base+os.sep+item
    #             out_dir_path = item.split(".")[0]
    #             out_dir_full_path=bs_cg_txt_base+os.sep+out_dir_path
    #             print "now handle %s" % out_dir_path
    #             out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
    #             pattern_type=TYPE_BS_CG_TXT
    #             bed_data_extract_to_methy(input_bed_path,out_bed_path,pattern_type)
    #             out_corr_dir=out_dir_full_path+os.sep+"correlation"
    #             # 如果文件夹不存在则应先创建
    #             if not os.path.exists(out_corr_dir):
    #                 os.makedirs(out_corr_dir)
    #             chr_no_list=range(1,23)
    #             for chr_i in chr_no_list:
    #                 bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
    #                 our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
    #                 calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)
    #
    wig_type_base="Homosapiens"+os.sep+"wig_type"
    wig_type_list=os.listdir(wig_type_base)

    for item in wig_type_list:
        if os.path.isfile(os.path.join(wig_type_base,item)):
            if item.split(".")[-1]=="wig":
                input_bed_path=wig_type_base+os.sep+item
                out_dir_path = item.split(".")[0]
                out_dir_full_path=wig_type_base+os.sep+out_dir_path
                if not os.path.exists(out_dir_full_path):
                    print "now handle %s" % out_dir_path
                    out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
                    variable_wig_data_extract_to_methy_normal(input_bed_path,out_bed_path)
                    out_corr_dir=out_dir_full_path+os.sep+"correlation"
                    # 如果文件夹不存在则应先创建
                    if not os.path.exists(out_corr_dir):
                        os.makedirs(out_corr_dir)
                    chr_no_list=range(1,23)
                    for chr_i in chr_no_list:
                        bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
                        our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
                        calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)
    #
    # encode_bed_type_base="Homosapiens"+os.sep+"encode_bed_type"
    # encode_bed_type_list=os.listdir(encode_bed_type_base)
    #
    # for item in encode_bed_type_list:
    #     if os.path.isfile(os.path.join(encode_bed_type_base,item)):
    #         if item.split(".")[-1]=="bed":
    #             input_bed_path=encode_bed_type_base+os.sep+item
    #             out_dir_path = item.split(".")[0]
    #             out_dir_full_path=encode_bed_type_base+os.sep+out_dir_path
    #             print "now handle %s" % out_dir_path
    #             out_bed_path=out_dir_full_path+os.sep+"splitted_bed"
    #             pattern_type=TYPE_ENCODE_BED
    #             bed_data_extract_to_methy(input_bed_path,out_bed_path,pattern_type)
    #             out_corr_dir=out_dir_full_path+os.sep+"correlation"
    #             # 如果文件夹不存在则应先创建
    #             if not os.path.exists(out_corr_dir):
    #                 os.makedirs(out_corr_dir)
    #             chr_no_list=range(1,23)
    #             for chr_i in chr_no_list:
    #                 bed_file = out_bed_path+os.sep+ "chr" + str(chr_i) + ".bed"
    #                 our_rd_path = out_corr_dir +os.sep+"chr"+ str(chr_i) + "_experiment_without.csv"
    #                 calc_correlation(chr_i, bed_file, our_rd_path, d_max, is_inter_with_other_cpg,ignore_d=ignore_d)