

rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathProfiles = "../input/PCAB_Ostalb_GD-pruned_MICA-per-tree"

# path where sample data is stored
pathSamples = "./samples"

########################################################



source("../libraries/mica-functions.R")

initMica(normalizePath("../libraries/"))

######################################
# read and store density profile data
######################################

yearMin = 1916
yearMax = 2004
rowsPerYear = 100
D = data.frame( rep(1:rowsPerYear,length(yearMin:yearMax)))
for (y in yearMin:yearMax) { D[(1:rowsPerYear)+(rowsPerYear*(y-yearMin)),2] = y}
colnames(D) = c("iInYear","year")

# get all tree files
treeFiles = list.files( pathProfiles , pattern=".*\\.csv")
# get available max density data from each file
for( tf in treeFiles ) {
  
  d = read.csv(paste(pathProfiles ,tf,sep="/"))
  
  Dcol = rep(NA, nrow(D))
  
  for (y in yearMin:yearMax) {
    if (sum(d$year == y) == rowsPerYear) {
      Dcol[ D$year == y ] = d$GD[ d$year == y ]
    }
  }
  
  D[,ncol(D)+1] = Dcol
  colnames(D)[ncol(D)] = substr( tf, 1, nchar(tf)-4)
  
} # for all treeFiles


######################################
#' distance between two profiles based on z-score slopes
#' @return average point-wise absolute slope differences of z-scores
getDistanceOfProfiles <- function( p1, p2 ) 
  ######################################
{
  if (length(p1)< 3) {
    stop(paste(p1))
  }
  s1 = getSlope( as.double(1:length(p1)), as.vector(scale(p1)) )
  s2 = getSlope( as.double(1:length(p2)), as.vector(scale(p2)) )
  # average point-wise absolute slope differences of z-scores
  return( mean( abs( s1$y - s2$y ) ) )
}

######################################
#' computes the distances for each putative start year of the sample within the chronology
#' @param chronoD the consensus profile information for each year
#' @param sampleD the profile information of the sample
#' @return the respective distance for each start year = sum(min profileDist per year)
getDistances <- function ( chronoD, sampleD ) 
  ######################################
{
  # initialize distances
  n = length(chronoD)/rowsPerYear
  d = rep(NA, n)
  nS = length(sampleD)/rowsPerYear
  for ( i in 1:(n-nS+1) )  {
    # overall distance = sum of minimal distance per year
    dist = 0;
    # for all years relative to current start year i
    for ( j in 0:(nS-1) ) {
      # find distance for this year
      jDist = getDistanceOfProfiles( p1=chronoD[((i-1+j)*rowsPerYear)+(1:rowsPerYear)]
                                     , p2=sampleD[(j*rowsPerYear)+(1:rowsPerYear)] )
      # add year's distance to overall distance
      dist = dist + jDist
    }
    # store overall distance
    d[i] = dist
  }
  return( d )
}


######################################
#' provides the row indices of the sample's profile within the overall chronology data
#' @param sample the sample encoding (file, start, end)
#' @return the profile data for the sample
getSampleRows <- function( sample ) {
  ######################################
  return( (1+ (sample$start - yearMin)*rowsPerYear):(rowsPerYear + (sample$end - yearMin)*rowsPerYear) )
}

######################################
#' extracts the sample's profile from the overall chronology data
#' @param chronoD the profile data per tree
#' @param sample the sample encoding (file, start, end)
#' @return the profile data for the sample
getSampleProfile <- function( chronoD, sample ) {
  ######################################
  return( chronoD[ getSampleRows( sample ), colnames(chronoD) %in% sample$file ] ) 
}


######################################
# get ranks in chronology for each sample set
######################################

for (l in c(10,15)){#,20)){#,5)) {
  
  
  # get all sample sets
  sampleFiles = list.files( pathSamples , pattern=paste("samples-len-",l,"-set-\\d+\\.csv",sep="") )
  
  allranks = c()
  
  # generate chronology for the set if not already present
  for (f in sampleFiles) {
    
    # read samples
    samples = as.data.frame(read.csv( paste(pathSamples,f,sep="/") ))
    
    # read samples' chronology
    chronoFile = paste(pathSamples,"/",substr(f,1,nchar(f)-4),".chronoConsensus.csv",sep="")
    if (!file.exists(chronoFile)) { next } # skip if no chronology available
    chrono = read.csv( chronoFile )
    
    # compute rank per sample
    sRank = rep(NA, nrow(samples))
    for ( s in 1:nrow(samples) ) {
      # get data of sample
      sD = getSampleProfile( D, samples[s,] )
      # compute dist for each shift in chronology
      cDist = getDistances( chrono$consGD, sD )
      # get rank for each shift
      cRank = rank(cDist)
      # store rank of sample
      sRank[s] = cRank[ (yearMin:yearMax) %in% samples$start[s]  ]
    } # for each sample
    
    # store ranks of samples
    allranks = cbind(allranks,sRank)
    
  } # for each sampleFile
  
  if (length(allranks)==0) { next } # skip if no data
  
  #sum(allranks==1)/sum(!is.na(allranks))
  sink("outfile.txt",append=TRUE)
  print(paste("consensus",
              "&",l,
              "&",
              round(mean(apply(allranks==1,2,sum)) / nrow(allranks)*100,digits=1),
              "&",
              round(mean(apply(allranks,2,median)),digits=1),
              "&",
              round(mean(apply(allranks,2,mean)),digits=1),
              "&",
              round(mean(apply(allranks,2,var),digits=1)),
              "\\"
  ))
  sink()
  #summary(allranks)
  
} # for sample length
