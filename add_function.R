args <- commandArgs(trailingOnly = TRUE)

# use shell input
snpp <- as.character(args[1])
table1<- as.character(args[2])
output_f<-as.character(args[3])

snp<-read.table(snpp, sep = ";", header = F)
table<-read.delim(table1)


#output folder of the new table and set it as wd
path <- output_f
setwd(path)

#print input and output folders
print(paste0("The sample is going to be processed:",as.character(args[2])))

##colnames of the processed gbk file
colnames(snp)<-c("gene_name","scaffold","Start","End","strand","Product")

##create intermediate table
start_row=c()
end_row=c()
Scaffolds=unique(snp$scaffold)
j=0
for(s in Scaffolds){
  j=j+1
  if(j%%10000==0){print(j/length(Scaffolds))}
  temp=which(snp$scaffold==s)
  x0=min(temp)
  xF=max(temp)
  start_row=c(start_row,x0)
  end_row=c(end_row,xF)
}
Int_table=data.frame(Scaffolds,start_row,end_row)

##new method
#vectors to become the new columns in the SNPs table
start_T=Sys.time()
gene=c()
product=c()
strand=c()
#dim(table)[1]
for (i in 1:dim(table)[1]){
  if(((i/1000)*100)%%1000==0){print(paste(i," out of ",dim(table)[1]," done! (",round(i/dim(table)[1]*100, digits = 1),"%)",sep=""))}
  
  s=as.character(table$scaffold[i])
  p=table$position[i]
  
  x=which(Int_table$Scaffolds==s)
  if(length(x)!=1){gene=c(gene, NA)
  product=c(product,NA)
  strand=c(strand,NA)}else{x0=Int_table$start_row[x]
  xF=Int_table$end_row[x]
  x=which(p>=snp$Start[x0:xF] & p<=snp$End[x0:xF])+x0-1
  if(length(x)!=1){gene=c(gene, NA)
  product=c(product,NA)
  strand=c(strand,NA)
  }else{gene=c(gene, as.character(snp$gene_name[x]))
  product=c(product,as.character(snp$Product[x]))
  strand=c(strand,as.character(snp$strand[x]))
  }
  }
  
}
end_T=Sys.time()
print(end_T-start_T)
##new table
new<-cbind(table, gene, product, strand)

name_table<-gsub(".tsv","", basename(table1))

write.csv(new, paste0(name_table,".csv"), sep=",", col.names = T, row.names = F)



