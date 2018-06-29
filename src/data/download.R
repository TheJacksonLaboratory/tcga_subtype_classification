
cohort.dat <- read.table("src/data/cohort.txt",sep="\t",as.is=T,header=T)

url.base1 <- "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/"
url.base2 <- "/20160128/gdac.broadinstitute.org_"
url.base3 <- ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz"

data_directory = "data/raw/"
for(i in 1:nrow(cohort.dat))
{
	cohort = cohort.dat[i,2]
	url <- paste(url.base1,cohort,url.base2,cohort,url.base3,sep="")
	dest.file <- paste(data_directory,cohort,".tar.gz",sep="")
	download.file(url,dest.file)
}


