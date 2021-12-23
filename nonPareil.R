###load package, if it is not present, it will install 
if(!require(Nonpareil)){install.packages("Nonpareil", repos = "http://cran.us.r-project.org")}

require(Nonpareil)

# get the input passed from the bash script: input of nonPaireil program and output the reports
args <- commandArgs(trailingOnly = TRUE)

# use shell input
path_input <- as.character(args[1])
path_output<- as.character(args[2])

#print input and output folders
print(paste0("The input directory was:", path_input ))
print(paste0("The output directory will be:", path_output))

##create the strings that will the /path/to/file.npo 
path <- file.path(path_input)
npo<-sort(list.files(path,pattern=".npo"))
pathfiles<-paste0(path,npo)

##set working directory from nonpareil report output foler from bash script
#empty list to store nonpareil reports
setwd(path_output)
pareilcurves<-list()

for (i in seq_along(pathfiles)){
  fil<-basename(pathfiles[i])
  file<-sub(".npo","", fil)
  file2<-sub("M_t-","", file)
  pdf(paste0(file2,".pdf"), height = 6, width = 7)
  curve<-Nonpareil.curve(pathfiles[i], star = 60, label = file2)
  pareilcurves[[i]]<- summary(curve)
  write.table(pareilcurves[i], paste0(file2,".txt"), sep = "\t", col.names = "F")
  dev.off()
}

print("nonPAreil curve Plots and reports done")


