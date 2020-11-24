import json
import csv
import os

with open('/home/raymond/Desktop/BlueSphere_Bio_Inc/MC38/MC38_BAM_somatic_mutation_calling/MC38_vcf_snpeff/keyword_filter_ANN_datafile/splice_region_variant.json') as f:
    snpeff_ann = json.load(f)
  
        
#clean the snpeff annotation file
snpeff_ann_dict = {}
for i in snpeff_ann.keys():
    snpeff_ann_dict[i.split(':')[0] + ':' + snpeff_ann[i][4]['ANN'][0].split('|')[6] + ':' +  i.split(':')[1]] = [snpeff_ann[i][0], snpeff_ann[i][1]]
print(list(snpeff_ann_dict.keys())[0], snpeff_ann_dict[list(snpeff_ann_dict.keys())[0]])
print(len(snpeff_ann_dict))
print(snpeff_ann_dict)


#check the feathers of the dataset(dictionary)
feather_count = {}
junk_residual = {}
for i in snpeff_ann_dict.keys():
    if len(snpeff_ann_dict[i][0]) == 1 and len(snpeff_ann_dict[i][1][1:-1]) == 1:
        feather_count[snpeff_ann_dict[i][0] + '->'+ snpeff_ann_dict[i][1]] = (len(snpeff_ann_dict[i][0]), len(snpeff_ann_dict[i][1][1:-1]))
    else:
        junk_residual[snpeff_ann_dict[i][0] + '->'+ snpeff_ann_dict[i][1]] = (len(snpeff_ann_dict[i][0]), len(snpeff_ann_dict[i][1][1:-1]))
print(feather_count)
print(len(feather_count))

#print(junk_residual)
#print(len(junk_residual))


#check if two snp positions are close in one transcript id 

#function: snp_close_check 
#format of input dataset: {'chromosome id: transcription id: snp position': ['ref, '[alt]']} 
#e.g., {'chr1:NM_001310513.1:36369105': ['C', '[A]']}
#the input dataset should follow this format otherwise the function won't work

def check_snp_close(input_dict):
    snp_close_dict = {}
    min_dis_dict = {}
    for k in input_dict.keys():
        if (k.split(':')[0], k.split(':')[1]) not in snp_close_dict.keys():
            snp_close_dict[(k.split(':')[0], k.split(':')[1])] = []
        snp_close_dict[(k.split(':')[0], k.split(':')[1])].append(k.split(':')[2])
        
    for j in snp_close_dict.keys():
            
        if len(snp_close_dict[j]) < 2:
            return 'No double snp in one transcription ID'
            
        if len(snp_close_dict[j]) == 2:
            snp_close_list = sorted(snp_close_dict[j])
            dis = int(snp_close_list[1]) - int(snp_close_list[0])
            if dis <= 12:
                return (k.split(':')[1], snp_close_list[1], snp_close_list[0])
            else:
                return 'There are two snps in one transcription ID but distance > 12', k.split(':')[0], k.split(':')[1]
            
        if len(snp_close_dict[j]) > 2:
            snp_close_list = sorted(snp_close_dict[j])
            for i in range(len(snp_close_list)-1):
                dis = int(snp_close_list[i+1]) - int(snp_close_list[i])
                if (snp_close_list[i+1], snp_close_list[i]) not in min_dis_dict:
                    min_dis_dict[(snp_close_list[i+1], snp_close_list[i])] = dis
                    if min(list(min_dis_dict.values())) <= 12:
                        min_dis_sort = sorted(list(min_dis_dict.values()))
                            
                    else:
                        return 'There are more than two snps in one transcription ID but distance > 12'
    return 'more than two snps in one transcription ID ', k.split(':')[0], k.split(':')[1], 'dis:' +  str(min_dis_sort[0])

#test:

#result1 = check_snp_close(snpeff_ann_dict)
#print(result1)

#check_dict2 = {'chr1:NM_001310513:15': ['C', '[A]'], 
#              'chr1:NM_001310513:20': ['T', '[A]']}
#result2 = check_snp_close(check_dict2)
#print(result2)

check_dict3 = {'chr1:NM_001310513:15': ['C', '[A]'], 
               'chr1:NM_001310513:20': ['T', '[A]'],
               'chr1:NM_001310513:22': ['G', '[C]'],
               'chr2:NM_001310514:10': ['C', '[A]'], 
               'chr2:NM_001310514:20': ['T', '[A]'],
               'chr2:NM_001310514:22': ['G', '[C]']
              
              }

result3 = check_snp_close(check_dict3)
print(result3)


file_path = '/home/raymond/Desktop/BlueSphere_Bio_Inc/MC38/MC38_BAM_somatic_mutation_calling/MC38_vcf_snpeff/keyword_filter_ANN_datafile/'

for i, j, file_name in os.walk(file_path):
    print(i, j, file_name)
    print(file_name)
    print(len(file_name))


result_list = []
for i in file_name:
    if i == 'splice_acceptor_variant.json'or i == 'splice_donor_variant.json'or  i == 'disruptive_inframe_insertion.json':
        del file_name[file_name.index(i)]
print(len(file_name))
            
for i in file_name:
    file = '/home/raymond/Desktop/BlueSphere_Bio_Inc/MC38/MC38_BAM_somatic_mutation_calling/MC38_vcf_snpeff/keyword_filter_ANN_datafile/{}'.format(i)
    with open (file) as f:
        snpeff_ann = json.load(f)
        
        snpeff_ann_dict = {}
        for j in snpeff_ann.keys():
            snpeff_ann_dict[j.split(':')[0] + ':' + snpeff_ann[j][4]['ANN'][0].split('|')[6] + ':' +  j.split(':')[1]] = [snpeff_ann[j][0], snpeff_ann[j][1]]
            result = check_snp_close(snpeff_ann_dict)
        result_list.append(result)
print(len(result_list)) 
print(result_list)

#check the features of result_list
for i in list(set(result_list)):
  print(i)
