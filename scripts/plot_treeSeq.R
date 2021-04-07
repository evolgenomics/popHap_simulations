args <- commandArgs(trailingOnly = TRUE)
tag<-args[1]

library(RColorBrewer)
colx<-brewer.pal(12,"Paired")
pdf(paste0("simTrees/PC062_sims.treeSeqs.",tag,".pdf"))
treeSeq<-read.table(paste0("simTrees/PC062_merged_Herat1603_3.45Mb.",tag,".treeSeq"), sep="\t", header=T)
plot(c(0,13504),c(1,968),xlab="SNP", ylab="Hap", main=paste0("Haplotype graph - ", tag),type="n");
cName<-colnames(treeSeq)
startStop<-lapply(cName, function(x) {strsplit(x,'.',fixed=T)[[1]][2:3]}) 
for (int in 1:length(startStop)) {
	for (i in 1:nrow(treeSeq)) {
		lines(as.numeric(startStop[[int]]), rep(i, 2), col=colx[(treeSeq[i,int]+1) %% 12]) 	
	}
}
dev.off()
