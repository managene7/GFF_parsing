#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':""}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print ("""
_____________________________________________________________________________

Usage;
python fasta_parsing_v2.0.py -in <file name> -out <file name> -db <file name> 

-help           show option list
-gff            input gff file name
-posi           position file name 
                name    start end (start and end are optional)
-out            output csv name (csv format)

_____________________________________________________________________________
""")
                quit()


def posi_parsing(posi_file):
    posi=open(posi_file,'r')
    posi_cont=posi.readlines()
    posi_dic={}
    for cont in posi_cont:
        cont=cont.strip()
        cont=cont.split()
        if cont[0] not in posi_dic:
            posi_dic[cont[0]]=[[int(cont[1]), int(cont[2])]]
        else:
            posi_dic[cont[0]].append([int(cont[1]), int(cont[2])])
    return posi_dic

def gff_parsing(gff_file, posi_dic,out_name):
    import csv
    out_csv=csv.writer(open(out_name,'w', newline=''))
    out_csv.writerow(['Query_name','Query_length','Query_from','Query_to','Hit_description','Hit_length','Hit_from','Hit_to','Score','e-value','Identity','Function'])
    gff=open(gff_file,'r')
    while 1:
        line=gff.readline()
        if line=="":
            print ("Completed!!")
            break
        line=line.strip()
        det=0
        if line[0] !="#":
            line=line.split('\t')
            if line[0] in posi_dic:
                name=line[0]
                posi_type=line[2]
                l_posi=int(line[3])
                r_posi=int(line[4])

                for posi in posi_dic[line[0]]:
                    if l_posi>=posi[0] and l_posi<=posi[1]:
                        det=1
                    if r_posi>=posi[0] and r_posi<=posi[1]:
                        det=1
                    if det==1:
                        if posi_type=="mRNA":
                            info=line[8].split(";")
                            for cont in info:
                                if "Name=" in cont:
                                    gene_name=cont[5:]
                                if "Note=" in cont:
                                    function=cont[6:-1]
                        det=0
                        if posi_type=="exon":
                            if l_posi>=posi[0] and l_posi<=posi[1]:
                                det=1
                            if r_posi>=posi[0] and r_posi<=posi[1]:
                                det=1
                            if det==1:
                                if l_posi<=posi[0]:
                                    cont_l_posi=1
                                else:
                                    cont_l_posi=l_posi-posi[0]+1
                                if r_posi>=posi[1]:
                                    cont_r_posi=posi[1]
                                else:
                                    cont_r_posi=r_posi-posi[0]+1
                                contents=[name,str(posi[1]-posi[0]+1), str(cont_l_posi), str(cont_r_posi), gene_name,str(cont_r_posi-cont_l_posi+1),'1', str(cont_r_posi-cont_l_posi+1),'100','0','1',function]
                                out_csv.writerow(contents)
                            det=0



def main():
    position=posi_parsing(option_dict['-posi'])
    print (position)
    if option_dict['-out']=='':
        option_dict['-out']=option_dict['-posi']+"_gff_parsed.csv"
    gff_parsing(option_dict['-gff'], position, option_dict['-out'])


if __name__=="__main__":
    main()