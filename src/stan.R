
library(rstan)
library(ggmcmc)
library(tidyverse)
library(rgl)
options(mc.cores = parallel::detectCores())

source("src/maphdi.R")


argv=commandArgs(T)
file1 =argv[1]
file2= argv[2]



##recipient 
#file1 <-"/home/ynanya/analysis/SNPs/tmp/GIFU-039-1calls2.txt"
name1 <-gsub("calls2.txt$", "", basename(file1))

data1 <-read_tsv(file1)

num1_A_allele <- unlist(data1[paste0(name1,"_A")])
num1_B_allele <- unlist(data1[paste0(name1,"_B")])


##donor
#file2 <-"/home/ynanya/analysis/SNPs/tmp/GIFU-039-1calls2.txt"
name2 <-gsub("calls2.txt$", "", basename(file2))

data2 <-read_tsv(file2)

num2_A_allele <- unlist(data2[paste0(name2,"_A")])
num2_B_allele <- unlist(data2[paste0(name2,"_B")])




##  determine recipient genotype 

freq1_A <- num1_A_allele /(num1_A_allele + num1_B_allele)
genotype1<- ifelse(freq1_A > 0.9, 1, 
                  ifelse(freq1_A < 0.6 & freq1_A > 0.4, 2,
                  ifelse(freq1_A<0.1, 3, 0)))
 #### genotype1 is a numerical vector 0:ND, 1:AA, 2:AB, 3:BB

valid = (!is.na(genotype1) & genotype1 >0 )
num_valid=sum(valid)

valid_num2_A_allele <- num2_A_allele[valid]
valid_num2_B_allele <- num2_B_allele[valid]
valid_genotype1 <-genotype1[valid]


valid2 = (valid_num2_A_allele + valid_num2_B_allele)>100
num_valid2=sum(valid2)

valid2_num2_A_allele = valid_num2_A_allele[valid2]
valid2_num2_B_allele = valid_num2_B_allele[valid2]
valid2_genotype1=valid_genotype1[valid2]




data <-list(N=num_valid2, valid_num2_A_allele=valid2_num2_A_allele, 
            total_num_allele= valid2_num2_A_allele + valid2_num2_B_allele,
            genotype=valid2_genotype1, theta=c(1,1,1))
fit <-stan(file='src/chimerism.stan', data=data, seed=1234)


ms <-rstan::extract(fit)


## draw 
pdf(paste0("result/",name1,"_",name2,".pdf"), width=10, height=16)
par(mfrow=c(3,2))
plot(unlist(num1_A_allele), unlist(num1_B_allele), xlab="A allele", ylab="B allele", main=name1)
plot(unlist(num2_A_allele), unlist(num2_B_allele), xlab="A allele", ylab="B allele", main=name2)

freq1_A <- as.numeric(unlist(num1_A_allele/(num1_A_allele + num1_B_allele)))
freq2_A <- as.numeric(unlist(num2_A_allele/(num2_A_allele + num2_B_allele)))

plot(c(1,0),c(0,1), col="white", main="A allele frequency", xlab="timing", ylab="A allele frequency")
for(s in 1:length(freq1_A)){
  lines(c(0,1), c(freq1_A[s], freq2_A[s]), type="l")
}


maphdi(fit, pars=c("chi"))-> chi_range
round(chi_range, 3) -> chi_range
hist(ms$chi[1001:2000], xlim=c(0,1), 
     main=(paste0("Histgram of donor chimerism\n", "chimerism:  ", chi_range["MAP"],":(", chi_range["Lower95"]," - ", chi_range["Upper95"], ")")),
     xlab="donor chimerism")

donor_genotype <- matrix(rep(0, 3*num_valid2), ncol=3) 
colnames(donor_genotype) <- c("AA", "AB", "BB")
for(n in 1:num_valid2){
  table(ms$donor_genotype[1001:2000,n]) -> tab
  donor_genotype[n,1] <- tab["1"]
  donor_genotype[n,2] <- tab["2"]
  donor_genotype[n,3] <- tab["3"]
  
}
barplot(t(donor_genotype), stack=T)


plot(c(-1,1),c(-0.2,2), col="white", main="Donor genotype estimation", 
     xlab="", ylab="", btn="n", xaxt="n", yaxt="n", axes=F)
polygon(c(-1,1,0), c(0,0,sqrt(3)))
text(-1, -0.07, "AA", cex=2)
text( 1, -0.07, "BB", cex=2)
text( 0, sqrt(3)+0.07, "AB", cex=2)

for(n in 1:num_valid2){
  table(ms$donor_genotype[1001:2000,n]) -> tab
  x =  0 + (tab["1"]*(-1) + tab["2"]*0 +tab["3"]*1) / 600
  y =  sqrt(3)/3.0 + (tab["1"]*(-1/sqrt(3)) + tab["2"]*(2/sqrt(3)) + tab["3"]*(-1/sqrt(3)) ) / 600 
  points(x,y)
  
}

par(mfrow=c(1,1))
dev.off()
#plot3d(donor_genotype[,1], donor_genotype[,2], donor_genotype[,3])




#out.file <- file("result/summary_result.txt", open = "a")
result_essence=paste(name1,name2,chi_range["MAP"], chi_range["Lower95"], chi_range["Upper95"], sep="\t")
#writeLines(result_essence, out.file)
#close(out.file)

write(result_essence, "result/summary_result.txt", append=T)

