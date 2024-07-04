#! /usr/bin/env Rscript


########################################################################
library(ape)


########################################################################
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = args[2]

if (length(args) < 2){ 
        print ("At least two arguments should be given! Exiting ......")
        q() 
}


########################################################################
MyTree <- read.tree(infile)

write.nexus(MyTree, file=outfile)


