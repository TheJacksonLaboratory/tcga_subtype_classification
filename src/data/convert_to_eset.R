
library(Biobase)

data_directory = "data/raw/"
zipped.files <- list.files(data_directory,pattern="gz",full=T)
for(i in 1:length(zipped.files))
{
	cmd <- paste("tar zxvf ",zipped.files[i],sep="")
	system(cmd)
}

gdac.dirs <- setdiff(list.dirs(),".")

for(i in 1:length(gdac.dirs))
{
	infile <- list.files(gdac.dirs[i],pattern="RSEM",full=T)
	rsem.file <- list.files(gdac.dirs[i],pattern="RSEM")
	type <- strsplit(rsem.file,".",fixed=T)[[1]][1]
	exp <- read.table(infile,sep="\t",header=F,skip=2,row.names=1)
	hdr <- read.table(infile,sep="\t",header=T,row.names=1,nrows=2)
	colnames(exp) <- colnames(hdr)
	entid <- sapply(rownames(exp),function(dat){ strsplit(as.character(dat),"|",fixed=T)[[1]][2]})
	symbol <- sapply(rownames(exp),function(dat){ strsplit(as.character(dat),"|",fixed=T)[[1]][1]})
	probe.anns <- data.frame(uiqid = rownames(exp),entid,symbol,stringsAsFactors=F)
	rownames(probe.anns) <- entid
	rownames(exp) <- entid
	f.dat <- new("AnnotatedDataFrame",probe.anns)
	cnames = gsub(".","-",substr(colnames(exp),1,15),fixed=T)
	sample.type = gsub(".","-",substr(colnames(exp),14,15),fixed=T)
	s.anns <- data.frame(sampleid = colnames(exp),sample.type,tumor.type = type,stringsAsFactors=FALSE)
	rownames(s.anns) <- cnames
	colnames(exp) <- cnames
	p.dat <- new("AnnotatedDataFrame",s.anns)
	eset = ExpressionSet(as.matrix(exp),phenoData = p.dat,featureData = f.dat)
	outfile <- paste(data_directory,type,"eset.rdata",sep="")
	save(eset,file = outfile)
}

for(i in 1:length(gdac.dirs))
{
	cmd <- paste("rm -rf ",gdac.dirs[i],sep="")
	system(cmd)
}
