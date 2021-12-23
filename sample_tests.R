args <- commandArgs(trailingOnly = TRUE)
# use shell input
path <- as.character(args[1])
setwd(path)

path1 <- file.path(path)
samples<-sort(list.files(path1,pattern="testt.samples"))

table<-read.table(samples)

##samples
cenas<-rep(c(1:1000), length.out=nrow(table), each= 2)
sample<-rep("Sample", times= length(cenas))
Samples<-paste0(sample,cenas)

#pairs
pair<-c("pair1","pair2")
pairs<-rep(pair, times= (length(cenas)/2))

#col bind to final table
table_final<-cbind(Samples,table,pairs)

##write table
write.table(table_final, "test.samples", 
            sep = "\t", 
            row.names = F, col.names = F,
            quote = F)


