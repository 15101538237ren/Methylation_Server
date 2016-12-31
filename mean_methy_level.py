# -*- coding:utf-8 -*-
import os,re

def calc_mean_methy_level(bedfile_dir,bedfile_path_list,out_file_path):
    out_file=open(out_file_path,"w")
    re_pattern = r'(\d+)\s([\d]+\.[\d]*)\s'
    for bed_file_path in bedfile_path_list:
        bed_file = open(bedfile_dir+os.sep+bed_file_path+".bed", 'r')
        line = bed_file.readline()
        cnt = 0
        sum_methy=0.0
        print "begin %s" % bed_file_path
        while line:
            match = re.search(re_pattern, line)
            if match:
                # pos = int(match.group(1))
                methy_level = float(match.group(2))
                sum_methy=sum_methy+methy_level
                cnt=cnt+1
            line = bed_file.readline()
        bed_file.close()
        print "end %s" % bed_file_path
        mean_methy_level=sum_methy/float(cnt)

        out_line=bed_file_path+","+str(mean_methy_level)+"\n"
        out_file.write(out_line)
    out_file.close()
    print "wrt succ!"
if __name__ == '__main__':
    bed_dir="Mouse_experiment_bed/chr1"
    path_list=["sperm","oocyte","2cell_paternal","2cell_maternal","4cell_paternal","4cell_maternal","E65_paternal","E65_maternal","E75_paternal","E75_maternal","E135M","E135F"]
    out_file_path="diff_period_mean_methy_level.csv"
    calc_mean_methy_level(bed_dir,path_list,out_file_path)