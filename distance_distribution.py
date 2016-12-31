# -*- coding:utf-8 -*-
import random,math,os,threading,datetime,re,shutil,collections,pickle,time
import numpy as np

#pickle 文件的生成见Promoter_CpG_Distribution.py 的calc_promoter_cpg_distribution函数
def calc_distance_distribution(dna_seq_pkl_dir,out_file_dir,chr_min=1,chr_max=20,strand_standard="+"):
    for chr_no in range(chr_min,chr_max):
        dna_seq_pkl_path=dna_seq_pkl_dir+os.sep+"chr"+str(chr_no)+".pkl"
        out_file_path=out_file_dir+os.sep+"chr"+str(chr_no)+".csv"
        out_file=open(out_file_path,"w")
        input_pickle = open(dna_seq_pkl_path, 'rb')
        print "loading %s!" % dna_seq_pkl_path
        sequence_array=pickle.load(input_pickle)
        print "%s loaded!" % dna_seq_pkl_path
        input_pickle.close()

        distance_freq_hash={}
        pos_list=[]
        print "joining sequences!"
        full_seq=("").join(sequence_array)
        print "join-seq finish!"
        if strand_standard=="+":
            iter_cg = re.finditer(r'CG',full_seq)
            for it in iter_cg:
                (first_pos,end_pos)=it.span()
                pos_list.append(first_pos)
        elif strand_standard=="-":
            iter_gc = re.finditer(r'GC',full_seq)
            for it in iter_gc:
                (first_pos,end_pos)=it.span()
                pos_list.append(end_pos)
        else:
            pos_list2=[]
            iter_gc = re.finditer(r'GC',full_seq)
            for it in iter_gc:
                (first_pos,end_pos)=it.span()
                pos_list2.append(end_pos)
            for l in range(len(pos_list2)-1):
                dist=pos_list2[l+1]-pos_list2[l]
                if dist not in distance_freq_hash.keys():
                    distance_freq_hash[dist]=1
                else:
                    distance_freq_hash[dist]=distance_freq_hash[dist]+1
            iter_cg = re.finditer(r'CG',full_seq)
            for it in iter_cg:
                (first_pos,end_pos)=it.span()
                pos_list.append(first_pos)

        for i in range(len(pos_list)-1):
            dist=pos_list[i+1]-pos_list[i]
            if dist not in distance_freq_hash.keys():
                distance_freq_hash[dist]=1
            else:
                distance_freq_hash[dist]=distance_freq_hash[dist]+1
        print "sorting distances!"
        sorted_hash= sorted(distance_freq_hash.iteritems(), key=lambda d:d[0])
        print "sorted!"
        str_list=[]
        for j in range(len(sorted_hash)):
            (dist,freq)=sorted_hash[j]
            str_list.append(str(dist)+","+str(freq))
        out_str=("\n").join(str_list)
        out_file.write(out_str)
        print "out of chr%d finished!" % chr_no
        out_file.close()
def tofel(input_file_path,out_file_path,out_hash_pkl_path=None,load=False):
    input_file = open(input_file_path,'r')
    line=input_file.readline()
    out_file = open(out_file_path, 'w')
    if load and out_hash_pkl_path!=None:
        pkl_file=open(out_hash_pkl_path,"rb")
        word_hash=pickle.load(pkl_file)
        pkl_file.close()
    else:
        word_hash={}
    pattern=r'\d+\.\s+([a-z]+)([^\n]+)'
    counter=0

    while line:
        match = re.search(pattern, line)
        if match:
            word=match.group(1)
            chinese=match.group(2)
            if word not in word_hash:
                counter=counter+1
                word_hash[word]=chinese
                line_to_wrt=str(counter)+"\t"+word+"\t"+chinese+"\n"
                out_file.write(line_to_wrt)
        line=input_file.readline()
    if out_hash_pkl_path!=None and load==False:
        pkl_file=open(out_hash_pkl_path,"wb")
        pickle.dump(word_hash,pkl_file,-1)
        pkl_file.close()
    out_file.close()
    input_file.close()
if __name__ == '__main__':
    tofel("/Users/Ren/Desktop/foot.txt","/Users/Ren/Desktop/foot_new.txt","/Users/Ren/Desktop/tofel.pkl",True)
    # species="mouse_cpg_distribution"
    # strand_standard="+"
    # chr_min=1
    # chr_max=20
    #
    # dna_pkl_path=species+os.sep+"dna"+os.sep+"pickle"
    # out_file_dir=species+os.sep+"distance_distribution"+os.sep+strand_standard
    # if not os.path.exists(out_file_dir):
    #     os.makedirs(out_file_dir)
    #
    # calc_distance_distribution(dna_pkl_path,out_file_dir,chr_min=chr_min,chr_max=chr_max,strand_standard=strand_standard)