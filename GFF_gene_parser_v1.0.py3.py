# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 17:41:10 2021

@author: minkj
"""
#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-len':50, '-fasta':""}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print ("""
_____________________________________________________________________________

Usage;

-help       show option list
-gff        gff file name
-posi       file name of seq position
-len        length of flanking regions to be parsed (kb unit, default is 10 kb)
-fasta      name of gene seq file (option)
_____________________________________________________________________________
""")
                quit()
                

def gff_to_dic(gff_file_name): #parsing gene info only from GFF3 file
    gff=open(gff_file_name,'r')
    out_dic={}
    print ("Convert GFF info into dic format..\n")
    for line in gff:
        if line[0] != "#":
            split_line=line.split()
            if split_line[2]=="gene":
                gene_mid_loc=int(int(split_line[3])+int(split_line[4])/2)
                #print (split_line)
                try:
                    gene_id=split_line[8].split(";")[0].split(":")[1]
                except:
                    try:
                        gene_id=split_line[8].split(";")[0].split("gene-")[1]
                    except:
                        gene_id=split_line[8].split(";")[0].split("=")[1]
                
                #if "description=" in split_line[8].split(";")[2]:
                #    gene_info=[split_line[8].split(";")[2].split("=")[1]]
                #else:
                #    gene_info=["-"]
                

                gene_left=split_line[3]
                gene_right=split_line[4]
                if split_line[0] not in out_dic:
                    out_dic[split_line[0]]={}
                    out_dic[split_line[0]][gene_mid_loc]=[gene_id, gene_left, gene_right]#+gene_info
                else:
                    out_dic[split_line[0]][gene_mid_loc]=[gene_id, gene_left, gene_right]#+gene_info
                #print ([gene_id, gene_left, gene_right]+gene_info)
    return out_dic

def gff_parsing(parsing_file, gff_dic_input): #parsing genes using position info
    gff_dic=gff_dic_input
    position=open(parsing_file,'r')
    position_parsed=position.readlines()
    out_dic={}
    print ("GFF info parsing..\n")
    for line in position_parsed:
        line=line.strip().split()

        chr=line[0]
        f_len=int(option_dict['-len'])*1000
        if int(line[1])>=f_len:
            l_posi=int(line[1])-f_len
        else:
            l_posi=0
        r_posi=int(line[1])+f_len
        out_key=line[0]+"_"+str(l_posi)+"-"+str(r_posi)
        out_dic[out_key]=[]
        
        if chr in gff_dic:
            gene_locs=gff_dic[chr].keys()
            for loc in gene_locs:
                loc=int(loc)
                if loc>=l_posi and loc<=r_posi:
                    out_dic[out_key].append(gff_dic[chr][loc])
        else:
            print ("Error in chromosome name to be parsed!")
            chr_list=list(gff_dic.keys())
            print ("chromosome names should be; ")
            print (chr_list)
            quit()
    return out_dic

def gene_seq_parsing(fasta):
    out_dic={}
    print ("Convert gene seq into dic format.. \n")
    gene_file=open(fasta,'r')
    init=0
    while 1:
        line=gene_file.readline().strip()
        if line=="":
            seq="".join(cont_dic)
            out_dic[name]=seq
            break
        
        else:
            if ">" in line:
                if init==0:
                    name=line
                    cont_dic=[]
                    init=1
                else:
                    seq="".join(cont_dic)
                    out_dic[name]=seq
                    
                    name=line
                    cont_dic=[]

            else:
                cont_dic.append(line)
    return out_dic
        

def main():
    import csv
    csv_out=csv.writer(open(option_dict['-posi']+"_parsed_gene_info.csv", 'w', newline=""))
    csv_out_list=[]
    csv_out.writerow(["Chr_l-loc_r-loc","Gene_ID", "gene_l-end", "gene_r-end"])
    gff_dic=gff_to_dic(option_dict['-gff'])
    parsed_dic=gff_parsing(option_dict['-posi'], gff_dic)

    for key, value in parsed_dic.items():
        if value==[]:
            cont=[key]+["No gene detected"]
            csv_out.writerow(cont)
            csv_out_list.append(cont)
        else:
            for sub_value in value:
                cont=[key]+sub_value
                csv_out.writerow(cont)
                csv_out_list.append(cont)
    if option_dict['-fasta']!="":
        print ("Parsing gene sequence..\n")
        gene_seq_dic=gene_seq_parsing(option_dict['-fasta'])
        gene_keys=gene_seq_dic.keys()
        
        gene_file=open(option_dict['-posi']+"_parsed_seqs.txt",'w')
        for row in csv_out_list:
            if row[0] != "Chr_l-loc_r-loc":
                gene_id=row[1]
                for gene_key in gene_keys:
                    if gene_id in gene_key:
                        gene_file.write(gene_key+"\n")
                        gene_file.write(gene_seq_dic[gene_key]+"\n")



if __name__=="__main__":
    main()



                
          