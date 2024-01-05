import os
import re

## For quality analysis, the software FastQC is required
## Run the command：
## python raw_fastqc.py

# Create a directory to save the quality analysis results, fastqc_result
os.system("mkdir fastqc_result")

# Only sequencing files with the suffix.fq in the current directory can be analyzed
os.system("ls *.fq > list_seq.txt")

file=open('list_seq.txt','r')
for i in file:
    fq=re.sub('\n','',i)
    # FastQC performed quality analysis of the sequencing files
    os.system("fastqc -o ./fastqc_result -t 5 %s" %(fq))
file.close()
os.system("rm list_seq.txt")

