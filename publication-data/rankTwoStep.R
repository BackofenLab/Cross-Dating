

rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathProfiles = "../input/PCAB_Ostalb_GD-pruned_MICA-per-tree"

# path where MICA alignment meta information is stored
pathMaxDens = "../input/PCAB_Ostalb_GD-max_densities"

# path where sample data is stored
pathSamples = "./samples"

rerankTopX = 20

########################################################


source("../libraries/mica-functions.R")

initMica(normalizePath("../libraries/"))



######################################
# read and store max density data
######################################

yearMin = 1916
yearMax = 2004
MD = data.frame(yearMin:yearMax)
colnames(MD) = c("year")

# get all tree files
treeFiles = list.files( pathMaxDens , pattern=".*\\.csv")
# get available max density data from each file
for( tf in treeFiles ) {
	
	d = read.csv(paste(pathMaxDens ,tf,sep=.Platform$file.sep))
	
	MD[ MD[,"year"] %in% d$year, ncol(MD)+1 ] = d$density
	colnames(MD)[ncol(MD)] = substr( tf, 1, nchar(tf)-4)
	
	
} # for all treeFiles


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
  
  d = read.csv(paste(pathProfiles ,tf,sep=.Platform$file.sep))
  
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
#' computes the distances for each putative start year of the sample within the chronology
#' @param chronoMD the max density chronology
#' @param sampleMD the sample's max densities
#' @return the distance for each start year
getDistancesMD <- function ( chronoMD, sampleMD ) {
	######################################
	# initialize distances
	n = length(chronoMD)
	d = rep(NA, n)
	nS = length(sampleMD)
	for ( i in 1:(n-nS+1) )  {
		d[i] = - cor( x=chronoMD[i:(i+nS-1)], y=sampleMD, method="pearson" )
	}
	return( d )
}


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
#' @param chronoD the profile information for each tree
#' @param sampleD the profile information of the sample
#' @param startYears either NA or a boolean vector of length nrow/rowsPerYear that marks whether the year is a valid start year
#' @return the respective distance for each start year = sum(min profileDist per year)
getDistances <- function ( chronoD, sampleD, startYears=NA ) 
######################################
{
  # initialize distances
  n = nrow(chronoD)/rowsPerYear
  d = rep(NA, n)
  nS = length(sampleD)/rowsPerYear
  for ( i in 1:(n-nS+1) )  {
    if ( !is.na(startYears) && !startYears[i] ) { next } # skip non-valid start years
    # overall distance = sum of minimal distance per year
    dist = 0;
    # for all years relative to current start year i
    for ( j in 0:(nS-1) ) {
      # get available profiles for year i+j (not NA in first pos of current year)
      jcols = !is.na(chronoD[ 1+ (i-1+j)*rowsPerYear, 3:ncol(chronoD)])
      jMinD = 999999; # init minimal distance for this year
      # find minimal distance for this year
      for ( c in 2+(1:length(jcols))[jcols] ) {
        jMinD = min( jMinD, getDistanceOfProfiles( p1=chronoD[((i-1+j)*rowsPerYear)+(1:rowsPerYear),c]
                                                   , p2=sampleD[(j*rowsPerYear)+(1:rowsPerYear)] ) )
      }
      # add year's minimal distance to overall distance
      dist = dist + jMinD
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
#' generates a chronology excluding the samples
#' @param chronoD the original chronology to prune
#' @param samples the samples to remove
#' @return the reduced chronology without sample data
getSampleChronology <- function( chronoD, samples ) {
######################################
  # copy chronology
  chrono = chronoD
  # set all samples' data to NA
  for (s in 1:nrow(samples)) {
    chrono[ getSampleRows( samples[s,] ), colnames(chrono) %in% samples$file[s] ] = NA
  }
  return(chrono)
}

######################################
# get ranks in chronology for each sample set
######################################

for (l in c(5)){#},10,15,20)) {


  # get all sample sets
  sampleFiles = list.files( pathSamples , pattern=paste("samples-len-",l,"-set-\\d+\\.csv",sep="") )
  
  allranks = c()
  
  # generate chronology for the set if not already present
  for (f in sampleFiles) {
    
    # read MD chrono data
    chronoFile = paste(substr(f,1,nchar(f)-4),"chronoMaxDens","csv",sep=".")
    if (!file.exists(paste(pathSamples,chronoFile,sep=.Platform$file.sep))) { stop(paste("chronofile",chronoFile,"not found")) }
    chronoMD = as.data.frame(read.csv( paste(pathSamples,chronoFile,sep=.Platform$file.sep), header=T ))
    
    # read samples
    samples = as.data.frame(read.csv( paste(pathSamples,f,sep=.Platform$file.sep) ))
    # generate sample chronology
    chrono = getSampleChronology( D, samples )
    
    # compute rank per sample
    sRank = rep(NA, nrow(samples))
    for ( s in 1:nrow(samples) ) {
      
      # (step 1) rank by max density
      
      # get data of sample
      sMD = MD[ MD$year %in% samples$start[s]:samples$end[s], colnames(MD) %in% samples$file[s] ]
      # compute dist for each shift in chronology
      cDist = getDistancesMD( chronoMD$maxDens, sMD )
      # get rank for each shift
      cRank = rank(cDist)
      
      # (step 2) rerank by density profiles
      
      # get top-ranked start years
      topRanked = cRank <= rerankTopX
      
      # get data of sample
      sD = getSampleProfile( D, samples[s,] )
      # compute dist for each shift in chronology
      cDist = getDistances( chrono, sD, startYears=topRanked )
      # get rank for each shift
      cRank[topRanked] = rank(cDist)[topRanked]
      # store rank of sample
      sRank[s] = cRank[ (yearMin:yearMax) %in% samples$start[s]  ]
    } # for each sample
    
    # store ranks of samples
    allranks = cbind(allranks,sRank)
    
  } # for each sampleFile


  #sum(allranks==1)/sum(!is.na(allranks))
  sink("outfile.txt",append=TRUE)
  print(paste("two step",
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
