
find.tumor.specific.features <- function(exp,sample.type, type="ACC")
{
	idx1 <- which(sample.type == type)
	idx2 <- which(sample.type != type)	
	
	t.test.fun <- function(x, idx1,idx2)
	{
		return(t.test(x[idx1],x[idx2])$p.value)
	}


	logFC.fun <- function(x, idx1,idx2)
	{
		return(mean(x[idx1]) -  mean(x[idx2]))
	}
	logFC <- apply(exp,1,logFC.fun,idx1,idx2)
	pvalue <- apply(exp,1,t.test.fun,idx1,idx2)
	res <- data.frame(pvalue,logFC)
	rownames(res) <- rownames(exp)
	res.sel <- res[res$pvalue < 0.001,]
	res.ord <- res.sel[order(res.sel$logFC,decreasing=T),]
	return(res.ord)
}


 
select_features_for_training <- function(learn.dat,num=25)
{
	cl.idx <- which(colnames(learn.dat) == "tumor.type")
	exp <- t(learn.dat[, -cl.idx])
	class.labels <- learn.dat[,cl.idx]
	types <- unique(class.labels)
	sel.features <- c()
	for(i in 1:length(types))
	{

		res <- find.tumor.specific.features(exp,class.labels, type=types[i])
		sel.features <- union(sel.features,rownames(res[1:num,]))
	}
	return(sel.features)
}



