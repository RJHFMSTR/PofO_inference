library(RColorBrewer)
COL = brewer.pal(8,"Set1")


D = read.table("test2", head=FALSE, stringsAsFactors=FALSE)
nD = ncol(D)

for (i in 1:11) {
	lb = paste("test_",i,".pdf", sep="")
	pdf(lb, 10, 5)
	par(mfrow=c(2,1))
	plot(0,0, type="l", xlim=c(0, nrow(D)), ylim=c(0,1), xlab="SNP index" , ylab="Copying prob.", main="First hap")
	points(D[, (i-1)*6+4], type="l", col=COL[1], lwd=2)
	points(D[, (i-1)*6+5], type="l", col=COL[2], lwd=2)
	points(D[, (i-1)*6+6], type="l", col="grey", lwd=2)
	legend("topright", legend=c("Father","Mother","Unrelated"), fill=c(COL[1], COL[2], "grey"), bg="white", cex=0.5)	
	plot(0,0, type="l", xlim=c(0, nrow(D)), ylim=c(0,1), xlab="SNP index" , ylab="Copying prob.", main="Second hap")
	points(D[, (i-1)*6+7], type="l", col=COL[1], lwd=2)
	points(D[, (i-1)*6+8], type="l", col=COL[2], lwd=2)
	points(D[, (i-1)*6+9], type="l", col="grey", lwd=2)
	dev.off()
}






