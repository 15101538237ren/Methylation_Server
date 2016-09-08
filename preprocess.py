# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle
import numpy as np
def bed_data_extract_to_methy(chr_no,in_file_path,outfile_path):
    #读入的文件路径
    raw_file = open(in_file_path,'r')

    #文件样例格式:chr1  3025349	3025350	0.6	3	2.染色体编号,CpG起始位点,CpG结束位点,甲基化水平,甲基化的reads数,未甲基化的reads数

    #写入的文件路径
    methy_file = open(outfile_path,'w')

    count=0
    find=0

    chr_no=str(chr_no)

    line=raw_file.readline()
    while line:
        count=count+1
        if count % 10000==0:
            print "%d lines of %s was processed" %(count,in_file_path)
        match = re.search(r'chr'+chr_no+r'\t(\d+)(\s)(\d+)(\s)([\d]+\.[\d]*)(\s)(\d+)(\s)(\d+)',line)
        if match and find==0:
            methy_pos=int(match.group(1))
            value=float(match.group(5))

            out_str=str(methy_pos)+"\t"+str(value)+"\n"

            methy_file.write(out_str)
            find=1
        elif match and find==1:
            methy_pos=int(match.group(1))
            value=float(match.group(5))
            out_str=str(methy_pos)+"\t"+str(value)+"\n"
            methy_file.write(out_str)
        #后边没有在match的CpG位点了
        elif find==1 and not match:
            break
        line=raw_file.readline()
    methy_file.close()
    raw_file.close()
    print "finish %s chr%s data processing!" %(in_file_path,chr_no)
def batch_bed_to_methy(input_file_path,chr_no_list,out_dir):
    #如果文件夹不存在则应先创建
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #遍历染色体列表
    for chr_no in chr_no_list:
        #每条染色体执行一遍提取bed数据操作输出到out_put_file_path
        out_put_file_path=out_dir+os.sep+"chr"+str(chr_no)+".bed"
        bed_data_extract_to_methy(str(chr_no),input_file_path,out_put_file_path)
if __name__ == '__main__':
    input_bed_path="GSM1386022_4cell_mc_CG_maternal_plus.bed"
    chr_no_list = [1]
    out_dir_path = "splitted_bed"
    batch_bed_to_methy(input_bed_path,chr_no_list, out_dir_path)