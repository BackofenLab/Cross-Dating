

source("../libraries/mica-functions.R")

initMica(normalizePath("../libraries/"))




# path where sample data is stored
pathSamples = "./samples"

# get all sample sets
micaFiles = list.files( pathSamples , pattern=paste("samples-len-\\d+-set-\\d+\\.chronoMICA.csv",sep="") )

rowsPerYear = 100

for (f in micaFiles) {

  # read data  
  d <- read.csv("samples/samples-len-10-set-1.chronoMICA.csv",header=T)
  
  # correct data
  d100=as.data.frame(matrix(NA,ncol=2, nrow=0))
  colnames(d100) = c("year","consGD")
  for (y in unique(d$year)) {
    last = nrow(d100)
    d100[last+(1:rowsPerYear),1] = rep(y,rowsPerYear)
    d100[last+(1:rowsPerYear),2] = interpolateCurve( as.double(1:sum(d$year==y)), d$GD[d$year==y], rowsPerYear )$y
  }
  
  # write data
  write.csv(d100, paste(pathSamples,"/", substr(f,1,nchar(f)-8),"Consensus.csv", sep=""),
            row.names = FALSE, quote = FALSE)
  
} # for all micaFiles



