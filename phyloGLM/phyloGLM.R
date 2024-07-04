library('phylolm')


#######################################################################
args <- commandArgs(trailingOnly = TRUE)
treefile = args[1]
datafile = args[2]
if (length(args) >= 3){
	isBin = TRUE
}else{
	isBin = FALSE
}


#######################################################################
tre = read.tree(file = treefile)
dat = read.table(file = datafile, header=T)

if (isBin){
	dat$predictor[dat$predictor >= 1] <- 1
}

#write.tree(tre)

fit = phyloglm(trait~predictor,phy=tre,data=dat,boot=100)
summary(fit)
coef(fit)
vcov(fit)
