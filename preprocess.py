# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle
import numpy as np
def bed_data_extract_to_methy(chr_no,in_file_path,outfile_path,standard=False):
    #读入的文件路径
    raw_file = open(in_file_path,'r')

    #文件样例格式:chr1  3025349	3025350	0.6	3	2.染色体编号,CpG起始位点,CpG结束位点,甲基化水平,甲基化的reads数,未甲基化的reads数

    #写入的文件路径

    if standard:
        out_plus_path=outfile_path+"_plus.bed"
        out_minus_path = outfile_path + "_minus.bed"
        methy_file1 = open(out_plus_path, 'w')
        methy_file2 = open(out_minus_path, 'w')
    else:
        methy_file1 = open(outfile_path,'w')

    count=0
    find=0

    chr_no=str(chr_no)

    line=raw_file.readline()
    while line:
        count=count+1
        if count % 100000==0:
            print "%d lines of %s was processed" %(count,in_file_path)
        if standard:
            match_p = re.search(r'chr' + chr_no + r'\t(\d+)\t(\d+)\t\.\t(\d+)\t\+\t(\d+)\t(\d+)\t(\d+),(\d+),(\d+)\t(\d+)\t(\d+)', line)
            match_m = re.search(r'chr' + chr_no + r'\t(\d+)\t(\d+)\t\.\t(\d+)\t-\t(\d+)\t(\d+)\t(\d+),(\d+),(\d+)\t(\d+)\t(\d+)', line)
            if match_p:
                methy_pos = int(match_p.group(1))
                value = float(match_p.group(10)) / 100.0
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
                if find == 0:
                    find = 1
            elif match_m:
                methy_pos = int(match_m.group(1))
                value = float(match_m.group(10)) / 100.0
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file2.write(out_str)
                if find == 0:
                    find = 1
            # 后边没有在match的CpG位点了
            elif find == 1 and not match_p and not match_m:
                break
        else:
            match_p = re.search(r'chr' + chr_no + r'\t(\d+)(\s)(\d+)(\s)([\d]+\.[\d]*)(\s)(\d+)(\s)(\d+)',line)
            if match_p:
                methy_pos = int(match_p.group(1))
                value = float(match_p.group(5))
                out_str = str(methy_pos) + "\t" + str(value) + "\n"
                methy_file1.write(out_str)
                if find == 0:
                    find = 1
            #后边没有在match的CpG位点了
            elif find==1 and not match_p:
                break
        line=raw_file.readline()
    methy_file1.close()
    if standard:
        methy_file2.close()
    raw_file.close()
    print "finish %s chr%s data processing!" %(in_file_path,chr_no)
def batch_bed_to_methy(input_file_path,chr_no_list,out_dir,standard=False):
    #如果文件夹不存在则应先创建
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #遍历染色体列表
    for chr_no in chr_no_list:
        #每条染色体执行一遍提取bed数据操作输出到out_put_file_path
        out_put_file_path=out_dir+os.sep+"chr"+str(chr_no)
        bed_data_extract_to_methy(str(chr_no),input_file_path,out_put_file_path,standard=standard)
if __name__ == '__main__':
    input_bed_path="GSM1386027_E135M_mc_CG_plus.bed"
    chr_no_list = range(1,2)
    out_dir_path = "E135M"
    batch_bed_to_methy(input_bed_path,chr_no_list, out_dir_path,standard=False)