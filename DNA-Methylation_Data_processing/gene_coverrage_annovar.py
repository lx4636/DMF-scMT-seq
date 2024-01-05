import pandas as pd
import numpy as np
import os
import datetime as dt
import argparse

import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sc', nargs='+',
                    default=[],
                    help='Single cell file num')
parser.add_argument('-b', '--bulk', nargs='+',
                    default='[]',
                    help='Bulk cell file num')
parser.add_argument('-n', '--name', type=str,
                    default='anno',
                    help='output file name')
parser.add_argument('-p', '--path', type=str,
                    default='./',
                    help='path')
                    
args = parser.parse_args()
singlecell = args.sc
bulkcell = args.bulk
filename = args.name
path = args.path
print('singlecell', singlecell)


#from final_file.methylkit_work.file_base.R_dir2.CpG_region_calc import calc_meth_bulk
#Specify the sample name for analysis DMF-{}.txt Single cell and bulk samples are separated
file_name_sc = singlecell
file_name_bulk = bulkcell
#Specify the gene region to view
anno_type = ['CpGisland','Promoters', 'Shores']
#Specify the information to analyze (the name of the output file)
file_name = filename
out_dir = path


date = str(dt.date.today())
def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print ("analyses directory created")
mkdir(out_dir)
################################

def methylation_list(df):
    #df=df.loc[df['coverage']>=1000]
    df['pos'] = df['chr']+'-'+df['base'].map(str)
    return df
def methylation_pos(df, df_cpg):
#     df['posi'] = df['chr']+'-'+df['start'].map(str)+'-'+df['end'].map(str)
    df_cpg_list = pd.merge(df, df_cpg,how='left', on='pos')
    df_cpg_list.fillna(0)
    #print(df_cpg_list)
    return df_cpg_list
def calc_meth(df):

    met = df[df['freqC']>=90]
    #df[df['Met']<3] =0
    met_num = len(met)
    total = len(df)
    if total == 0:
        return 0
    methyl_level = (round((met_num/total)*100, 2))
    return methyl_level
def calc_meth_stat_bulk(df):
    df['MetC'] = df['coverage']*df['freqC']/100
    met_num = sum(df['MetC'])
    #df[df['Met']<3] =0
    total = sum(df['coverage'])
    if total == 0:
        return 0
    methyl_level = (round((met_num/total)*100, 2))
    return methyl_level
def Coverrate(file, CpG_site, Type):
    test_89 = pd.read_csv(file, sep='\t')
    test_89 = methylation_list(test_89)
    test_89_anno = methylation_pos(test_89, CpG_site)
    a = test_89_anno[test_89_anno['type']==Type]
    a = a[a['coverage']>=3]
    a['MetC'] = a['coverage']*a['freqC']/100
    return calc_meth(a), calc_meth_stat_bulk(a)
def methylation_level_calculation(df):
    df_intergenic = df[df['GeneType']=='intergenic']
    df_exonic = df[(df['GeneType'] == 'exonic')|(df['GeneType'] == 'UTR5')|(df['GeneType']=='UTR3')|(df['GeneType']=='UTR5;UTR3')]
    df_intronic = df[(df['GeneType'] == 'intronic')]
    df_splicing = df[(df['GeneType'] == 'splicing')]
    df_updown = df[(df['GeneType'] == 'upstream')|(df['GeneType'] == 'downstream')|(df['GeneType'] =='upstream;downstream')]
    intergenic_rate = calc_meth_stat_bulk(df_intergenic)
    exonic_rate = calc_meth_stat_bulk(df_exonic)
    intronic_rate =calc_meth_stat_bulk(df_intronic)
    splicing_rate = calc_meth_stat_bulk(df_splicing)
    updown_rate = calc_meth_stat_bulk(df_updown)
    return [intergenic_rate,exonic_rate,intronic_rate,splicing_rate,updown_rate]
def methylation_level_calculation_SC(df):
    df_intergenic = df[df['GeneType']=='intergenic']
    df_exonic = df[(df['GeneType'] == 'exonic')|(df['GeneType'] == 'UTR5')|(df['GeneType']=='UTR3')|(df['GeneType']=='UTR5;UTR3')]
    df_intronic = df[(df['GeneType'] == 'intronic')]
    df_splicing = df[(df['GeneType'] == 'splicing')]
    df_updown = df[(df['GeneType'] == 'upstream')|(df['GeneType'] == 'downstream')|(df['GeneType'] =='upstream;downstream')]
    intergenic_rate = calc_meth(df_intergenic)
    exonic_rate = calc_meth(df_exonic)
    intronic_rate =calc_meth(df_intronic)
    splicing_rate = calc_meth(df_splicing)
    updown_rate = calc_meth(df_updown)
    return [intergenic_rate,exonic_rate,intronic_rate,splicing_rate,updown_rate]
def vcf_input_generate(df_merge):
    ID = ['.' for i in range(len(df_merge['chr']))]
    score = [100 for i in range(len(df_merge['chr']))]
    filte = ['PASS' for i in range(len(df_merge['chr']))]
    INFO = ['DP=4' for i in range(len(df_merge['chr']))]
    vcf_input = pd.DataFrame({'#CHROM':df_merge['chr'], 'POS':df_merge['base'], 'ID':ID, 'REF':ID, 'ALT':ID, 'QUAL':score, 'FILTER':filte, 'INFO':INFO})
    return vcf_input

result_df = pd.DataFrame()
result_df['gene'] = ['Intergenic_rate','Exonic_rate','Intronic_rate','Splicing_rate','Updown_rate']
print('begin')

if 'CpGisland' in anno_type:
    test = pd.read_csv('/home/songjia/bigdisk/xx/reference_file/hg19_cpg.bed', sep='\t', names = ['chr','start', 'end', 'CpG'])
    chr_n = list(test['chr'])
    start = list(test['start'])
    end = list(test['end'])
    pos = []
    for i in range(len(chr_n)):
        pos.extend([chr_n[i]+'-'+str(j) for j in range(start[i], end[i]+1)])
    CpGisland_site = pd.DataFrame({'pos': pos})
    CpGisland_site['type'] = ['CpGisland' for i in range(len(CpGisland_site))]
if 'Shores' in anno_type:
    shore = pd.read_csv('/home/songjia/bigdisk/xx/final_file/methylkit_work/file_base/R_dir2/shores.txt', sep=',',index_col=0)
    chr_n = list(shore['seqnames'])
    start = list(shore['start'])
    end = list(shore['end'])
    pos = []
    for i in range(len(chr_n)):
        pos.extend([chr_n[i]+'-'+str(j) for j in range(start[i], end[i]+1)])
    Shores_site = pd.DataFrame({'pos': pos})
    Shores_site['type'] = ['Shores' for i in range(len(Shores_site))]
if 'Promoters' in anno_type:
    promoter = pd.read_csv('/home/songjia/bigdisk/xx/final_file/methylkit_work/file_base/R_dir2/promoter.txt', sep=',',index_col=0)
    chr_n = list(promoter['seqnames'])
    start = list(promoter['start'])
    end = list(promoter['end'])
    pos = []
    for i in range(len(chr_n)):
        pos.extend([chr_n[i]+'-'+str(j) for j in range(start[i], end[i]+1)])
    Promoters_site = pd.DataFrame({'pos': pos})
    Promoters_site['type'] = ['Promoters' for i in range(len(Promoters_site))]

if file_name_sc != '[]': 
    for i in file_name_sc:
        tmp = pd.read_csv(i, sep='\t')
        vcf_tem = vcf_input_generate(tmp)
        vcf_tem.to_csv(i.split('.')[0]+'.vcf', index=False, header=True, sep='\t')
        os.system('perl /home/songjia/bigdisk/xx/software/annovar/convert2annovar.pl --format vcf4 {} > {}'.format(i.split('.')[0]+'.vcf',i.split('.')[0]+'.avinput'))
        os.system('perl /home/songjia/bigdisk/xx/software/annovar/table_annovar.pl {} /home/songjia/bigdisk/xx/software/annovar/humandb -otherinfo --build hg19 -out {} -protocol refGene,cytoBand -operation g,r -remove -nastring \".\"'.format(i.split('.')[0]+'.avinput',i.split('.')[0]))
        print('annovar done!!!')
        print('Read the annotated file')
        ann_tmp = pd.read_csv(i.split('.')[0]+'.hg19_multianno.txt', sep='\t')
        tmp['GeneType'] = ann_tmp['Func.refGene']
        result = methylation_level_calculation_SC(tmp)
        result_df[str(i)] = result
if file_name_bulk != '[]':
    for i in file_name_bulk:
        tmp = pd.read_csv(i, sep='\t')
        vcf_tem = vcf_input_generate(tmp)
        vcf_tem.to_csv(i.split('.')[0]+'.vcf', index=False, header=True, sep='\t')
        os.system('perl /home/songjia/bigdisk/xx/software/annovar/convert2annovar.pl --format vcf4 {} > {}'.format(i.split('.')[0]+'.vcf',i.split('.')[0]+'.avinput'))
        os.system('perl /home/songjia/bigdisk/xx/software/annovar/table_annovar.pl {} /home/songjia/bigdisk/xx/software/annovar/humandb -otherinfo --build hg19 -out {} -protocol refGene,cytoBand -operation g,r -remove -nastring \".\"'.format(i.split('.')[0]+'.avinput',i.split('.')[0]))
        print('annovar done!!!')
        print('Read the annotated file')
        ann_tmp = pd.read_csv(i.split('.')[0]+'.hg19_multianno.txt', sep='\t')
        tmp['GeneType'] = ann_tmp['Func.refGene']
        result = methylation_level_calculation(tmp)
        result_df[str(i)] = result
result_df.to_csv(out_dir+'anno_'+file_name+date+'result_slice1.csv', sep=',')
#result_df = pd.read_csv('./anno_2021-11-09result_slice1.csv', sep=',')
if anno_type:
    result_df_2 = pd.DataFrame()
    result_df_2['gene'] = ['CpGisland','Promoters', 'Shores']
    if file_name_sc != '[]': 
        for j in file_name_sc:
            res_lis = []
            for i in anno_type:
                res,_ = Coverrate(str(j), locals()[i+'_site'], i)
                #print(j,'done',j, res)
                res_lis.append(res)
                print(i, 'done')
            print(res_lis)
            result_df_2[j] = res_lis
    if file_name_bulk != '[]':
        for j in file_name_bulk:
            res_lis = []
            for i in anno_type:
                _,res = Coverrate(str(j), locals()[i+'_site'], i)
                #print(j,'done',j, res)
                res_lis.append(res)
            print(res_lis)
            result_df_2[j] = res_lis
    result_df_2.to_csv('anno_'+date+'result_slice2.csv', sep=',')
    result_df = pd.merge(result_df, result_df_2, how='outer')    
result_df.to_csv(out_dir+'anno_test_'+file_name+date+'result.csv', sep=',',index=False)
d= pd.read_csv(out_dir+'anno_test_'+file_name+date+'result.csv', sep=',',index_col='gene')
d=d.T
d.to_csv(out_dir+'anno_test_'+file_name+date+'result.csv',sep=',')
d.to_csv(out_dir+'anno_noIndex_'+file_name+date+'result.csv',sep=',', index=False)

