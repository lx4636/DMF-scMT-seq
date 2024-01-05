# run quality control with trim_galore
file1=(*_1.fq.gz)
file2=(*_2.fq.gz)
for ((i=0;i<=${#file1[@]};i++)); do
trim_galore --quality 20 --phred33 --stringency 3 --gzip --length 36 --rrbs --paired --cores 8 --trim1 --output_dir ./cleandata ${file1[i]} ${file2[i]}
done
###########################
echo "quality control done"
# run lamda DNA alignment with bismark
file3=(./cleandata/*val_1.fq.gz)
file4=(./cleandata/*val_2.fq.gz)
for ((i=0;i<=${#file3[@]};i++)); do 
bismark /home/songjia/bigdisk/xx/reference/reference_bowtie1/ bowtie2 -quiet -o ./cleandata/lamdaRna_alignment/ --temp_dir ./cleandata/lamdaRna_alignment/tmp /home/songjia/bigdisk/xx/reference/reference_bowtie1/ -1 ${file3[i]} -2 ${file4[i]}
done

echo "lamda DNA alignment done"
# compute the conversion rate
for file in ./cleandata/lamdaRna_alignment/*.bam
do
    echo $file 
    samtools view -h $file > ${file/.bam/.sam}
done

for file in ./cleandata/lamdaRna_alignment/*.sam
do
    echo $file
    perl /home/songjia/bigdisk/xx/Xiamen-University/47-56/lamdaRna_alignment/MethylExtractBSCR.pl seqFile=/home/songjia/bigdisk/xx/reference/NC_001416.fa inFile=$file flagW=99,147 flagC=83,163 > ${file/.sam/.transrate} 
done

echo "trans_rate compute done"
# run  DNA alignment with bismark
file1=(./cleandata/*val_1.fq.gz)
file2=(./cleandata/*val_2.fq.gz)
for ((i=0;i<=${#file1[@]};i++)); do
bismark /home/songjia/bigdisk/xx/reference/hg19 bowtie2 -quiet -o ./cleandata/alignment_result_bowtie2/ --multicore 8 --temp_dir ./cleandata/alignment_result_bowtie2/tmp /home/songjia/bigdisk/xx/reference/hg19 -1 ${file1[i]} -2 ${file2[i]}
done

echo "DNA alignment done"

# run BAM sort with picard 
cd ./cleandata/alignment_result_bowtie2


echo "bam sort done"
# path generation
if [ ! -d ./sort_bam ]; then
  mkdir ./sort_bam
fi

for i in *.bam
do
    echo $i
    samtools sort -o ${i/.bam/_sorted.bam} -@ 4 $i
done
