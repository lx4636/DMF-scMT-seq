import pandas as pd 
import numpy as np 
import os
import subprocess
import argparse
import glob
import datetime
import re

##########Reading in parameters#############
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--state', type=str,
                    default='start',
                    help='start medium end')
parser.add_argument('-p', '--path', type=str,
                    default='./',
                    help='path')
parser.add_argument('-t', '--target', type=str,
                    default='./',
                    help='sorted bam')
parser.add_argument('-T', '--type', type=str,
                    default='RRBS',
                    help='WGBS RRBS')
                
args = parser.parse_args()
state = args.state

########Path Settings##########
path_default = args.path
path_sorted_bam = args.target

start_time = datetime.datetime.now()
time_str = start_time.strftime('%Y_%m_%d_%H_%M')

mode = args.type
if mode == 'RRBS':
    path_t = '/cleandata/alignment_result_bowtie2/'
else:
    path_t = '/cleandata/alignment_result_bowtie2/deduplicate'
path = os.getcwd()


#######Deduplication-alignment-sort###########
def control():
    if mode == 'RRBS':
        os.system('sh /home/songjia/bigdisk/xx/scripts/process.sh')
    else:
        os.system('sh /home/songjia/bigdisk/xx/scripts/WGBS-process.sh')
    
############Follow-up analysis############
def methyl_call():  
    if not os.path.exists('./tmp/'):
        os.makedirs('./tmp/')
    sorted_bam = []
    for sortbam in sorted(glob.glob('*sorted.bam')):
        print ("\n current file is: " + sortbam)
        sorted_bam.append(sortbam.split('/')[-1])
    bam_df = pd.DataFrame({'file': sorted_bam})
    a_path = './tmp/file_list'
    bam_df.to_csv(a_path, header=False, index=False)
    print('a_path', a_path)
    os.system("/home/songjia/software/anaconda/anaconda3/envs/r-4.0/bin/Rscript /home/songjia/newdisk/lx/scripts/DNAmethylation/Methylkit_call.R -i {} -o {} > Methylkit_result".format(a_path,'./'))
    
    CpG_file = []
    for infile in sorted(glob.glob('*_CpG.txt')):
        print ("\n current file is: " + infile)
        CpG_file.append(infile)
    os.system('python /home/songjia/newdisk/lx/scripts/DNAmethylation/gene_coverrage_annovar.py -s {} -p {}'.format(' '.join(CpG_file), './'))
    os.system('rm -r ./tmp')

##########Get coverage and conversion rates##############
def get_cover():
    def CpG_coverage(df_path):
        tmp = pd.read_csv(df_path, sep='\t')
        return len(tmp)
    cover_dict = {}
    for infile in sorted(glob.glob('*_CpG.txt')):
        print ("\n current file is: " + infile)
        cover_dict.update({infile: CpG_coverage(infile)})
    tmp_df = pd.DataFrame({'Sample': list(cover_dict.keys()), 'Cover': list(cover_dict.values())})
    tmp_df.to_csv('./'+time_str+'cover_result.csv')
    
    
def get_transrate():
    def get_rate(file):
        file_t = open(file, 'r+')
        file_txt = file_t.read()
        rate = float(re.findall(r'conversion rate: (.*?)\n', file_txt)[0])
        return rate
    transrate_dict = {}
    for infile in sorted(glob.glob('*.transrate')):
        print ("\n current file is: " + infile)
        transrate_dict.update({infile: get_rate(infile)})
    print(transrate_dict)
    tmp_df = pd.DataFrame({'Sample': list(transrate_dict.keys()), 'trans': list(transrate_dict.values())})
    tmp_df.to_csv('./'+time_str+'transrate_result.csv')   


  
if state == 'start':
    control()
    os.chdir(path+path_t)
    os.mkdir('{}'.format(time_str+'_report'))
    os.system('bismark2report --dir {}'.format(time_str+'_report'))
    methyl_call()
    get_cover()
    os.chdir('../lamdaRna_alignment')
    
    get_transrate()
if state == 'medium':
    methyl_call()
    get_cover()
    os.chdir('../lamdaRna_alignment')
    get_transrate()
if state == 'end':
    get_cover()
    os.chdir('../lamdaRna_alignment')
    get_transrate()

