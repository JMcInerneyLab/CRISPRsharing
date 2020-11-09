import os
import re
import glob
import errno


def nonBlankLine(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

#CRISPRStudio gff rewrite

d_filename = {} #create set for all sample file .fna name and its NC_** name

with open('filename_accession.txt') as f:
    for line in f:
        locus = line.split('\t')[0]
        acc = line.split('\t')[1]
        d_filename[locus] = acc


PATH = '/Users/Sera/Documents/NewCRISPR/Result_Analysis/Re_collection2/Gain_Loss/Sample*/Cluster*.txt'
filelist = glob.glob(PATH)

#copy fasta file into seperate directory
for i_a in filelist:
    print(i_a)
    dirc = i_a.split('/')[8][6:]
    d_cas_array = {} #dictionary for locus and its array number from CRISPRCasfinder
    d_det_array = {} #dictionary for locus and its array number from CRISPRDetect
    det_file = set() #create a set of .fna file who have CRISPRDetect results
    cas_file = set() #create a set of .gff file who have CRISPRCasfinder results
    renumber_array_det = {}
    with open(i_a) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                locus = '_'.join(line[1:].split('_')[:-2])
                if line.endswith('CAS'):
                    array_num = line.split('_')[-2]
                    if locus in d_cas_array.keys():
                        d_cas_array[locus].add(array_num)
                    else:
                        stock = set()
                        stock.add(array_num)
                        d_cas_array[locus] = stock
                    cas_file.add(locus)
                if line.endswith('DET'):
                    array_num = line.split('_')[-2]
                    if locus in d_det_array.keys():
                        d_det_array[locus].add(array_num)
                    else:
                        stock = set()
                        stock.add(array_num)
                        d_det_array[locus] = stock
                    det_file.add(d_filename[locus])
                    
    all_array_number_check = 0
    for j in d_cas_array.values():
        all_array_number_check += len(j)
    for j in d_det_array.values():
        all_array_number_check += len(j)
    print(all_array_number_check)
    
    #Filter arrays from CRISPRDetect results    
    PATH = '/Users/Sera/Documents/NewCRISPR/Result_Analysis/Re_collection2/CRISPR_GFF/det_gff_results/GCF*.fna.detect.txt.gff'

    filelist2S = glob.glob(PATH)
    filelist2= list()

    for i in det_file:
        for j in filelist2S:
            if i in j:
                filelist2.append(j)
                
    with open('Sample{}/CRISPRArrayResult{}.gff'.format(dirc,dirc),'w') as f_out:
        for i in filelist2:
            try:
                with open(i) as f:
                    for line in f:
                        locus = line.split('\t')[0].split('.')[0]
                        sarray = line.split('\t')[8].split(';')
                        renumber_array_det = {}
                        num = 1
                        if locus in d_det_array.keys():
                            for array_number_dic in sorted(d_det_array[locus]):
                                renumber_array_det[array_number_dic] = num
                                num+=1
                            for it in sarray:
                                if re.match('ID=',it):
                                    array_num = find_between(it,"ID=CRISPR","_")
                                    if array_num in d_det_array[locus]:
                                        f_out.write(line.replace('ID=CRISPR{}'.format(array_num),'ID=CRISPR{}'.format(renumber_array_det[array_num])))                   
            except IOError as exc:
                if exc.errno != errno.EISDIR:
                    raise

    #CRISPRCasFinder gff result collection                        
    with open('Sample{}/CRISPRArrayResult{}.gff'.format(dirc,dirc),'a') as f_out:
        for lo in cas_file:
            with open('/Users/Sera/Documents/NewCRISPR/Result_Analysis/Re_collection2/CRISPR_GFF/cas_gff_results/{}.gff'.format(lo)) as f:
                for line in nonBlankLine(f):
                    if line.startswith('N'):
                        Array = line.split('\t')
                        accession = Array[0]
                        software = Array[1]
                        atype = Array[2]
                        start_p = int(Array[3])
                        end_p = int(Array[4])
                        length = abs(start_p-end_p) +2
                        direction = Array[6]
                        renumber_array_cas={}
                        if renumber_array_det == {}:
                            num = 0
                        else:
                            num = sorted(renumber_array_det.values())[-1]
                        for array_number_dic in sorted(d_cas_array[accession.split('.')[0]]):
                            num+=1
                            renumber_array_cas[array_number_dic] = num 
                        if direction == '.':
                            direction = '+'
                        sarray = Array[8].split(';')
                        if atype == 'CRISPR':
                            c = 1
                            for i in sarray:
                                if re.match('DR=',i):
                                    sequence = i[3:]
                                if re.match('ID=',i):
                                    array_num = i.split('_')[-1]
                                    if array_num in d_cas_array[accession.split('.')[0]]:
                                        sArray = 'ID=CRISPR{}_{}_{};Note={};;Ontology_term=CRISPR'.format(renumber_array_cas[array_num],start_p,end_p+1,sequence)
                                        f_out.write('{}\tCRISPRCasfinder\trepeat_region\t{}\t{}\t{}\t{}\t.\t{}\n'.format(accession,start_p,end_p,length,direction,sArray))
                                        for line in f:
                                            Array = line.split('\t')
                                            accession = Array[0]
                                            software = Array[1]
                                            atype = Array[2]
                                            start = int(Array[3])
                                            end = int(Array[4])
                                            length = abs(start-end) +1
                                            direction = Array[6]
                                            if direction == '.':
                                                direction = '+'
                                            sarray = Array[8].split(';')
                                            for it in sarray:
                                                if re.match('sequence=',it):
                                                    sequence = it[9:]
                                            if atype == 'CRISPRdr':
                                                sArray = 'ID=CRISPR{}_REPEAT{}_{}_{};Name=CRISPR{}_REPEAT{}_{}_{};Parent=CRISPR{}_{}_{};Note={};;Ontology_term=CRISPR'.format(renumber_array_cas[array_num],c,start,end+1,renumber_array_cas[array_num],c,start,end+1,renumber_array_cas[array_num],start_p,end_p+1,sequence)
                                                f_out.write('{}\tCRISPRCasfinder\tdirect_repeat\t{}\t{}\t{}\t{}\t.\t{}\n'.format(accession,start,end,length,direction,sArray))
                                                
                                            if atype == 'CRISPRspacer':
                                                sArray = 'ID=CRISPR{}_SPACER{}_{}_{};Name=CRISPR{}_SPACER{}_{}_{};Parent=CRISPR{}_{}_{};Note={};;Ontology_term=CRISPR'.format(renumber_array_cas[array_num],c,start,end+1,renumber_array_cas[array_num],c,start,end+1,renumber_array_cas[array_num],start_p,end_p+1,sequence)
                                                f_out.write('{}\tCRISPRCasfinder\tbinding_site\t{}\t{}\t{}\t{}\t.\t{}\n'.format(accession,start,end,length,direction,sArray))
                                                c +=1
                                            if atype == 'RightFLANK':
                                                break

