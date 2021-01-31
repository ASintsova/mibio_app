# Title     : TODO
# Objective : TODO
# Created by: ansintsova
# Created on: 23.01.21

install.packages(c('BiocManager','Hmisc'),
                 dependencies='Depends',
                repo = "http://cran.us.r-project.org")
BiocManager::install(c('DESeq2','apeglm'))



library(DESeq2)

args <- commandArgs(trailingOnly=TRUE)

sdf_file <-args[1]

edf_file <- args[2]
experiment <- unlist(strsplit(basename(sdf_file), '_'))[[1]]

sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)
calculate_fitness <- function(sdf, edf){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=~day)
    dds <- DESeq2::DESeq(dds)
    fitness <- list()
    for(day in levels(dds$day)[-1]){
        fitness[[day]] = as.data.frame(DESeq2::results(dds,contrast=c("day",day,"d0")))
    }

    lapply(levels(dds$day)[-1], function(x) write.table(fitness[[x]], paste("../data/tmp_", experiment, '_', x, ".txt", sep="")))
}

calculate_fitness(sdf, edf)

