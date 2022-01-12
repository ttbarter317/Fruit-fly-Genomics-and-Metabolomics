# CMH tests on 10 Selection regime SNP table
# The permuation tests are run in parallel so be sure to adjust the number of cores used
library(foreach)
library(doParallel)
registerDoParallel(cores=60) # This should be set to the machine capabiliities 

# In the next statement put the correct file name in the read.table statement.
snp.data <- read.table("meta_cov20_snptable.txt", header=FALSE, col.names=c("chr", "pos", "major_nuc", "minor_nuc", "ac1_mjc", "ac1_mc", "ac2_mjc", "ac2_mc", "ac3_mjc", "ac3_mc", "ac4_mjc", "ac4_mc", "ac5_mjc", "ac5_mc","ao1_mjc", "ao1_mc", "ao2_mjc", "ao2_mc", "ao3_mjc", "ao3_mc", "ao4_mjc", "ao4_mc", "ao5_mjc","ao5_mc","co1_mjc","co1_mc", "co2_mjc", "co2_mc", "co3_mjc", "co3_mc", "co4_mjc", "co4_mc", "co5_mjc","co5_mc","nco1_mjc","nco1_mc", "nco2_mjc", "nco2_mc", "nco3_mjc", "nco3_mc", "nco4_mjc", "nco4_mc", "nco5_mjc","nco5_mc","scoa1_mjc","scoa1_mc","scoa2_mjc","scoa2_mc", "scoa3_mjc", "scoa3_mc", "scoa4_mjc", "scoa4_mc", "scoa5_mjc", "scoa5_mc","scob1_mjc", "scob1_mc", "scob2_mjc", "scob2_mc","scob3_mjc", "scob3_mc", "scob4_mjc","scob4_mc","scob5_mjc","scob5_mc"))
n<- length(snp.data[,1])

# Next we construct at each SNP a 10 x 2 contigency table over the five replicates
# and then apply the CMH test and save the p-values in "p.values".
p.values<- sapply(1:n,function(i){
  x<- snp.data[i,]
  snp_i_table<- array(c(x$ac1_mjc,x$co1_mjc,x$scoa1_mjc,x$ac1_mc,x$co1_mc,x$scoa1_mc,x$ac2_mjc,x$co2_mjc,x$scoa2_mjc,x$ac2_mc,x$co2_mc,x$scoa2_mc,x$ac3_mjc,x$co3_mjc,x$scoa3_mjc,x$ac3_mc,x$co3_mc,x$scoa3_mc,x$ac4_mjc,x$co4_mjc,x$scoa4_mjc,x$ac4_mc,x$co4_mc,x$scoa4_mc,x$ac5_mjc,x$co5_mjc,x$scoa5_mjc,x$ac5_mc,x$co5_mc,x$scoa5_mc,x$ao1_mjc,x$nco1_mjc,x$scob1_mjc,x$ao1_mc,x$nco1_mc,x$scob1_mc,x$ao2_mjc,x$nco2_mjc,x$scob2_mjc,x$ao2_mc,x$nco2_mc,x$scob2_mc,x$ao3_mjc,x$nco3_mjc,x$scob3_mjc,x$ao3_mc,x$nco3_mc,x$scob3_mc,x$ao4_mjc,x$nco4_mjc,x$scob4_mjc,x$ao4_mc,x$nco4_mc,x$scob4_mc,x$ao5_mjc,x$nco5_mjc,x$scob5_mjc,x$ao5_mc,x$nco5_mc,x$scob5_mc),dim=c(3,2,10))
  man.test<-mantelhaen.test(snp_i_table)
  return(c(man.test$p.value,man.test$statistic))
})

# Now align the p-values with the chromosome label and nucleotide position. The p.table$test.stat has a chisq
# distribution with 9 degrees of freedom (e.g. number of selection regimes -1)
# p.table can also be used for plotting the test statistic as a funciton of the genetic map
p.table<- data.frame(chr=snp.data$chr,pos=snp.data$pos,p.values=p.values[1,],test.stat=p.values[2,])

# Next we will calculate the false discovery rate (FDR) using the plug in method
# outlined in "Elements of Statistical Learning", pg. 683-690
# The plug in method requires comparing the test statistics and choosing cut points ahead of time
# to determine their FDR.
# We first compute the frequency of observations that exceed each cut point in the observations.
cut.points<- c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600)
Robs<- sapply(1:length(cut.points), function(i) length(p.table[p.table[,4]>cut.points[i],1]))

# Now do permutations and determine the frequency of observations greater than the cut points
n.perm<- 60
col.num<- (1:30)*2 +3 # creates of sequence from 5 to 103 skipping every other number
sim.cut<- foreach(i=1:n.perm,.errorhandling = c('pass')) %dopar%{
  col.sample<- sample(col.num,30)
  col.sample2<- c(1,2,3,4)
  for (j in 1:30) col.sample2<- c(col.sample2,col.sample[j],col.sample[j]+1)
  perm.data<- snp.data[,col.sample2]
  names(perm.data)<-c("chr", "pos", "major_nuc", "minor_nuc", "ac1_mjc", "ac1_mc", "ac2_mjc", "ac2_mc", "ac3_mjc", "ac3_mc", "ac4_mjc", "ac4_mc", "ac5_mjc", "ac5_mc","ao1_mjc", "ao1_mc", "ao2_mjc", "ao2_mc", "ao3_mjc", "ao3_mc", "ao4_mjc", "ao4_mc", "ao5_mjc","ao5_mc","co1_mjc","co1_mc", "co2_mjc", "co2_mc", "co3_mjc", "co3_mc", "co4_mjc", "co4_mc", "co5_mjc","co5_mc","nco1_mjc","nco1_mc", "nco2_mjc", "nco2_mc", "nco3_mjc", "nco3_mc", "nco4_mjc", "nco4_mc", "nco5_mjc","nco5_mc","scoa1_mjc","scoa1_mc","scoa2_mjc","scoa2_mc", "scoa3_mjc", "scoa3_mc", "scoa4_mjc", "scoa4_mc", "scoa5_mjc", "scoa5_mc","scob1_mjc", "scob1_mc", "scob2_mjc", "scob2_mc","scob3_mjc", "scob3_mc", "scob4_mjc","scob4_mc","scob5_mjc","scob5_mc")
  sim.tests<- sapply(1:n,function(j){
    x<- perm.data[j,]
    snp_i_table<- array(c(x$ac1_mjc,x$co1_mjc,x$scoa1_mjc,x$ac1_mc,x$co1_mc,x$scoa1_mc,x$ac2_mjc,x$co2_mjc,x$scoa2_mjc,x$ac2_mc,x$co2_mc,x$scoa2_mc,x$ac3_mjc,x$co3_mjc,x$scoa3_mjc,x$ac3_mc,x$co3_mc,x$scoa3_mc,x$ac4_mjc,x$co4_mjc,x$scoa4_mjc,x$ac4_mc,x$co4_mc,x$scoa4_mc,x$ac5_mjc,x$co5_mjc,x$scoa5_mjc,x$ac5_mc,x$co5_mc,x$scoa5_mc,x$ao1_mjc,x$nco1_mjc,x$scob1_mjc,x$ao1_mc,x$nco1_mc,x$scob1_mc,x$ao2_mjc,x$nco2_mjc,x$scob2_mjc,x$ao2_mc,x$nco2_mc,x$scob2_mc,x$ao3_mjc,x$nco3_mjc,x$scob3_mjc,x$ao3_mc,x$nco3_mc,x$scob3_mc,x$ao4_mjc,x$nco4_mjc,x$scob4_mjc,x$ao4_mc,x$nco4_mc,x$scob4_mc,x$ao5_mjc,x$nco5_mjc,x$scob5_mjc,x$ao5_mc,x$nco5_mc,x$scob5_mc),dim=c(3,2,10))
    man.test<-mantelhaen.test(snp_i_table)
    return(man.test$statistic)
  })
  V<- sapply(1:length(cut.points), function(j) length(sim.tests[sim.tests>cut.points[j]]))
  V
}

# Now we need to unlist sim.cut and make it a matrix
sim.cut.mtx<- NULL
for (i in 1:n.perm) sim.cut.mtx<- cbind(sim.cut.mtx,unlist(sim.cut[[i]]))

# Calculate the expected value of V
EV<- sapply(1:length(cut.points), function(i) sum(sim.cut.mtx[i,])/n.perm)
FDR<- EV/Robs #FDR has the false discovery rate for the CMH tests for each cut point in cut.points

p.done<-p.table[p.table[,4]>=500,]
p.done1<-p.done[p.done[,1]=="2L" | p.done[,1]=="2R" | p.done[,1]=="3L" | p.done[,1]=="3R" | p.done[,1]=="X",]
done<-merge(snp.data,p.done1,by.x=c("chr","pos"))
write.table(done,"meta_20cov_sigSNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
