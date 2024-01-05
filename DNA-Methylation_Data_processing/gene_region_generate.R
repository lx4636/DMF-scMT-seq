library(methylKit)
library(genomation)
file.list=list('DMF-59_CpG.txt','DMF-92_CpG.txt',
               'DMF-129_CpG.txt', 'DMF-130_CpG.txt')
myobj=methRead(file.list,
           sample.id=list('59','92','129','130'),
           assembly="hg19",
           treatment=c(1,1,0,0),
           context="CpG",
           mincov = 0
           )
meth=unite(myobj, destrand=FALSE)
file_name <- c('59','92','129','130')
gene.obj=readTranscriptFeatures("hg19.bed")
cpg.obj=readFeatureFlank("hg19_cpg.bed",
                           feature.flank.name=c("CpGi","shores"))
promoters=regionCounts(myobj,gene.obj$promoters)
exons=regionCounts(myobj,gene.obj$exons)
introns=regionCounts(myobj,gene.obj$introns)
Tss=regionCounts(myobj,gene.obj$TSSes)
Cpgi=regionCounts(myobj,cpg.obj$CpGi)
shores=regionCounts(myobj,cpg.obj$shores)
for (i in c(1,2,3,4)){
    write.csv(promoters[i], paste(file_name[i],'promoter.txt',sep='_'))
    write.csv(exons[i], paste(file_name[i],'exons.txt',sep='_'))
    write.csv(introns[i], paste(file_name[i],'introns.txt',sep='_'))
    write.csv(Tss[i], paste(file_name[i],'Tss.txt',sep='_'))
    write.csv(Cpgi[i], paste(file_name[i],'Cpgi.txt',sep='_'))
    write.csv(shores[i], paste(file_name[i],'shores.txt',sep='_'))
} 