
cat("input parameters:\n")
print(arg)

int_limit = as.numeric(arg[3])
dat<-read.table(arg[1],stringsAsFactors=F,sep="\t",comment.char="",head=T)
pdf(arg[2],width=6.5, height=3.5)
par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
for(i in 1:dim(dat)[1]){
	query = dat[i,]
	mz<-as.numeric(strsplit(split=";",x=as.character(query["mz"]))[[1]])
	int<-as.numeric(strsplit(split=";",x=as.character(query['intensity']))[[1]])
	int.max<-max(int,na.rm=T)
	int<-int/int.max*100
	m.label<-strsplit(split=";",x=as.character(query['m_label']))[[1]]
	m.mz<-as.numeric(strsplit(split=";",x=as.character(query['m_mz']))[[1]])
	m.int<-as.numeric(strsplit(split=";",x=as.character(query['m_intensity']))[[1]])
	m.int<-m.int/int.max*100

	## filter low intensity peak
	high_index <- m.int > int_limit
	m.label <- m.label[high_index]
	m.mz <- m.mz[high_index]
	m.int <- m.int[high_index]
	
	m.label=gsub(pattern="#|&",replacement="",x=m.label)
	plot(mz,int,type="h",col="gray",ylim=c(0,120),yaxs="i",cex.lab=1.01,font.lab=2,main=paste(query['Query'],query['pepSeq'],sep=" "),cex.main=0.65,xlab="MZ",ylab="Intensity(%)")
	if(length(m.label)>=1){
		text(m.mz,m.int,labels=paste(m.label,sprintf("%.4f",m.mz),sep=" "),cex=0.6,adj=c(-0.1,0.5),srt=90,col=ifelse(grepl(pattern="y",x=m.label),"red","blue"))
		lines(m.mz,m.int,type="h",lwd=1.1,col=ifelse(grepl(pattern="y",x=m.label),"red","blue"))
	}
	
}
dev.off()
