
# variable alert
ref=/home/songjia/reference/hg19/hg19_latest.fa


file1=(*_1.fastq.gz)
file2=(*_2.fastq.gz)
if [ ! -d no_fliter_cleandata ]; then
    mkdir no_fliter_cleandata
fi


for ((i=0;i<=${#file1[@]};i++)); do
    fileN=${file1[i]}
    trim_galore --quality 20 --phred33 --stringency 3 --gzip --length 36 --paired --cores 4 --trim-n --output_dir ./no_fliter_cleandata/${fileN/_1.fastq.gz} $fileN ${fileN/_1.fastq.gz}_2.fastq.gz
    cd ./no_fliter_cleandata/${fileN/_1.fastq.gz}
    file=${fileN/_1.fastq.gz}
    # if [ ! -d baw_align ]; then
    #     mkdir baw_align
    # f    # cd baw_align

    # if [ ! -d ${file1/_1.fq.gz} ]; then
    #     mkdir ${file1/_1.fq.gz}
    # fi
    bwa mem -t 5 -M $ref ${file}_1_val_1.fq.gz ${file}_2_val_2.fq.gz > ${file}.sam  

    # picard sort sam
    picard SortSam I=${file}.sam O=${file}_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR TMP_DIR=./tmp

    rm -r ./tmp

    if [ ! -d ${file}_snp ]; then
        mkdir ${file}_snp
    fi
    
    #mkdir ${file/.bam}_snp 
    
    cd ${file}_snp
    
    #bam files were marked duplicates with picard MarkDuplicates to avoid PCR duplicate reads affecting subsequent call SNPS.
    picard MarkDuplicates I=../${file}_sorted.bam O=${file}_marked.bam METRICS_FILE=KPGP-00001_L1.metrics
    #samtools mpileup
    samtools mpileup -go samtools.bcf -f $ref -t DP -t SP ${file}_marked.bam 
    #bcftools 
    bcftools call -vmO z -o bcftools_raw.vcf.gz samtools.bcf
    #Filter QUAL less than 10, DP values less than 5, and sites near the INDEL
    #bcftools filter -O v -o bcftools_filter.vcf -s LOWQUAL -e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . bcftools_raw.vcf.gz
    # Extracting SNP loci
    bcftools view -v snps bcftools_raw.vcf.gz > ${file}_bcftools_snp_raw.vcf
    # bcftools analyses the loci
    bcftools stats ${file}_bcftools_snp_raw.vcf >  ${file}_view.stats
    if [ ! -d ${file}_plot]; then
        mkdir ${file}_plot
    fi
    /home/songjia/software/bcftools-1.8/misc/plot-vcfstats ${file}_view.stats -p ./${file}_plot
    cd ../../..
done
