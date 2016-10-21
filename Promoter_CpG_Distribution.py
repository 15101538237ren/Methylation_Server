# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle,time
import numpy as np

def extract_all_TSS_pos_from_gtf(input_gtf_path,pickle_path,upstream=4000,downstream=2000,strand_standard="+",chr_min=1,chr_max=20,pattern=r'(\d+)\s+protein_coding\s+exon\s+(\d+)\s+(\d+)\s+[^\s]+\s+([-|+])\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+exon_number\s+\"1\"[^\n]+'):
    #读入的文件路径
    annotation_file = open(input_gtf_path,'r')
    line=annotation_file.readline()
    out_pickle = open(pickle_path, 'wb')

    count=0
    #up记录正链上游的启动子的截止位置
    promoter_up=[]
    #down记录正链下游启动子的截止位置
    promoter_down=[]

    for i in range(chr_max):
        promoter_up.append([])
        promoter_down.append([])

    str_pattern=pattern

    while line:
        count=count+1
        if count % 100000==0:
            print "%d lines of %s was processed" %(count,input_gtf_path)
        match_p = re.search(str_pattern, line)
        if match_p:
            chr_no=int(match_p.group(1))
            if chr_no <= chr_max and chr_no>= chr_min:
                strand=match_p.group(4)
                if strand==strand_standard:
                    if strand_standard=="+":
                        tss=int(match_p.group(2))
                        up=tss-upstream

                        down=tss+downstream
                    else:
                        tss=int(match_p.group(3))
                        up=tss-downstream
                        down=tss+upstream
                    promoter_up[chr_no].append(up)
                    promoter_down[chr_no].append(down)
        line=annotation_file.readline()
    annotation_file.close()


    pickle.dump(promoter_up,out_pickle,-1)
    pickle.dump(promoter_down,out_pickle,-1)
    print "%s dump completed!" % pickle_path
    out_pickle.close()

def get_row_and_col_by_index(index):
    row=index/60

    col=index%60
    return row,col
def calc_promoter_cpg_distribution(pickle_path,dna_seq_pkl_dir,out_file_path,sequence_dir,dna_file_pre,seq_storage=True,chr_min=1,chr_max=20,strand_standard="+"):
    input_pickle = open(pickle_path, 'rb')
    promoter_up=pickle.load(input_pickle)
    promoter_down=pickle.load(input_pickle)
    print "%s read completed!" % pickle_path
    input_pickle.close()
    out_file=open(out_file_path,"w")

    for chr_no in range(chr_min,chr_max):
        if not os.path.exists(dna_seq_pkl_dir):
            os.makedirs(dna_seq_pkl_dir)
        dna_seq_pkl_path=dna_seq_pkl_dir+os.sep+"chr"+str(chr_no)+".pkl"

        if seq_storage:
            out_pickle = open(dna_seq_pkl_path, 'wb')
            sequence_array=[]
            dna_file=open(sequence_dir+os.sep+dna_file_pre+str(chr_no)+".fa","r")
            line_seq=dna_file.readline()
            match = re.search(r'chromosome:([^:]+):(\d*):(\d*):(\d*)',line_seq)
            line_counts=int(math.ceil(int(match.group(4))/60.0+1.0))
            percent_100=line_counts/100
            line_no=1
            while line_no < line_counts+1:
                if (line_no % percent_100== 0):
                    percent = line_no/percent_100
                    print "%d percent for reading dna file of %d" % (percent,chr_no)
                line_seq = dna_file.readline()
                sequence_array.append(line_seq[0:60])
                line_no=line_no+1

            pickle.dump(sequence_array,out_pickle,-1)
            print "%s dump completed!" % dna_seq_pkl_path
            dna_file.close()
            out_pickle.close()
        else:
            input_pickle = open(dna_seq_pkl_path, 'rb')
            print "loading %s!" % dna_seq_pkl_path
            sequence_array=pickle.load(input_pickle)
            print "%s loaded!" % dna_seq_pkl_path
            input_pickle.close()

        promoter_up_arr=promoter_up[chr_no]
        promoter_down_arr=promoter_down[chr_no]

        promoter_len=promoter_down_arr[0]-promoter_up_arr[0]

        freq_of_cpg=[]
        for j in range(promoter_len):
            freq_of_cpg.append(0)
        first=True
        for i in range(len(promoter_up_arr)):
            up=promoter_up_arr[i]
            row_up,col_up=get_row_and_col_by_index(up)
            down=promoter_down_arr[i]
            row_down,col_down=get_row_and_col_by_index(down)

            #get promoter sequence
            if row_down < len(sequence_array):
                first_line_str=sequence_array[row_up][col_up:60]
                interval_line_str=""
                row_now=row_up+1
                while row_now < row_down:
                    interval_line_str=interval_line_str+sequence_array[row_now]
                    row_now=row_now+1
                last_line_str=sequence_array[row_down][0:col_down]
                full_promoter=first_line_str+interval_line_str+last_line_str
                if first:
                    iter_cg = re.finditer(r'CG',full_promoter)
                    for it in iter_cg:
                        (first_pos,end_pos)=it.span()
                        freq_of_cpg[first_pos]=freq_of_cpg[first_pos]+1
        freq_of_cpg_str=[]

        for k in range(promoter_len):
            freq_of_cpg_str.append(str(freq_of_cpg[k]))

        out_line_str=(",").join(freq_of_cpg_str)
        out_file.write(out_line_str+"\n")

        print "finish chr%d writiing" % chr_no
    out_file.close()
if __name__ == '__main__':
    species="zebrafish_cpg_distribution"
    strand_standard="-"
    chr_min=1
    chr_max=26

    seq_storage=False #DNA序列是否第一次读取,需要缓存成pkl

    pickle_path=species+os.sep+"promoter_tss_"+strand_standard+".pkl"
    dna_pkl_path=species+os.sep+"dna"+os.sep+"pickle"
    sequence_dir=species+os.sep+"dna"+os.sep+"sequences"

    human_dna_file_pre="Homo_sapiens.GRCh38.dna.chromosome."
    mouse_dna_file_pre="Mus_musculus.NCBIM37.67.dna.chromosome."
    zebrafish_dna_file_pre="Danio_rerio.Zv9.dna.chromosome."

    human_gtf_path=sequence_dir+os.sep+"Homo_sapiens.GRCh38.86.chr.gtf"
    mouse_gtf_path=sequence_dir+os.sep+"Mus_musculus.NCBIM37.67.gtf"
    zebrafish_gtf_path=sequence_dir+os.sep+"Danio_rerio.Zv9.79.gtf"

    human_pattern=r'(\d+)\s+[^\s]+\s+exon\s+(\d+)\s+(\d+)\s+[^\s]+\s+([-|+])\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+exon_number\s+\"1\"[^\n]+'
    mouse_pattern=r'(\d+)\s+protein_coding\s+exon\s+(\d+)\s+(\d+)\s+[^\s]+\s+([-|+])\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+exon_number\s+\"1\"[^\n]+'
    zebrafish_pattern=human_pattern


    out_file_path=species+os.sep+"chr"+str(chr_min)+"_"+str(chr_max-1)+strand_standard+".csv"



    # extract_all_TSS_pos_from_gtf(zebrafish_gtf_path,pickle_path,pattern=zebrafish_pattern,strand_standard=strand_standard,chr_min=chr_min,chr_max=chr_max)

    calc_promoter_cpg_distribution(pickle_path,dna_pkl_path,out_file_path,sequence_dir,zebrafish_dna_file_pre,seq_storage=seq_storage,chr_min=chr_min,chr_max=chr_max,strand_standard=strand_standard)
