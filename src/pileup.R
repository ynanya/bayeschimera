setwd("/home/ynanya/analysis/SNPs/chimerism/")


argv=commandArgs(T)
samplename =argv[1]
samplepath =argv[2]
work_dir =argv[3]
SNPs_path =argv[4]

	
dfTarget<-data.frame(specimen=as.character(samplename), bam=as.character(samplepath) )
path_tmp_file_Target<-paste0(work_dir,"/", samplename, ".Targetbam")
write.table(dfTarget, path_tmp_file_Target, col.names=TRUE, sep="\t", quote=FALSE,  row.names=FALSE)

php_command <- paste("php src/call_SNP.php", samplename, SNPs_path, work_dir, " 2", sep=" ")
print(php_command)
system ( php_command )

file.remove(path_tmp_file_Target)

