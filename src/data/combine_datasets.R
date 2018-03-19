library(Biobase)

data_directory = "../../data/raw/"

infiles <- list.files(pattern="rdata",data_directory,full=T)

load(infiles[1])
combined.eset <- eset
for(i in 2:length(infiles))
{
        load(infiles[i])
        combined.eset <- combine(combined.eset,eset)
	print(paste("i = ",i))
}

outfile <- paste(data_directory,"combined_data.rdata",sep="")
save(combined.eset, file = outfile)

