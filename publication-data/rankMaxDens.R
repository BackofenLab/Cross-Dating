

rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathMaxDens = "../input/PCAB_Ostalb_GD-max_densities"

# path where sample data is stored
pathSamples = "./samples"

########################################################


# implementations
# if(!exists("ANALYZER_IMPORTED")) source("../maths/Analyzer.R");
# if(!exists("CLUSTERING_IMPORTED")) source("../maths/Clustering.R");
# if(!exists("CURVE_MINER_IMPORTED")) source("../maths/CurvesMiner.R");
# if(!exists("DATA_ANALYSIS_IMPORTED")) source("../experiments/DataAnalysis.R");
# if(!exists("DEFAULTS_IMPORTED")) source("../system/Defaults.R");
# if(!exists("ENCODER_IMPORTED")) source("../visuals/Encoder.R");
# if(!exists("EXPERIMENT_0_IMPORTED")) source("../experiments/Experiment0.R");
# if(!exists("EXPERIMENT_1_IMPORTED")) source("../experiments/Experiment1.R");
# if(!exists("EXPERIMENT_2_IMPORTED")) source("../experiments/Experiment2.R");
# if(!exists("EXPERIMENT_3_IMPORTED")) source("../experiments/Experiment3.R");
# if(!exists("EXPERIMENT_4_IMPORTED")) source("../experiments/Experiment4.R");
# if(!exists("EXPERIMENT_5_IMPORTED")) source("../experiments/Experiment5.R");
# if(!exists("INTERFACE_IMPORTED")) source("../system/Interface.R");
# if(!exists("INTERPRETER_IMPORTED")) source("../system/Interpreter.R");
# if(!exists("LOADER_IMPORTED")) source("../system/Loader.R");
# if(!exists("MATH_IMPORTED")) source("../maths/Math.R");
# if(!exists("PLOTTER_IMPORTED")) source("../visuals/Plotter.R");
# if(!exists("PRESENTATION_IMPORTED")) source("../experiments/Presentation.R");
# if(!exists("STORER_IMPORTED")) source("../system/Storer.R");


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
  
  d = read.csv(paste(pathMaxDens ,tf,sep="/"))
  
  MD[ MD[,"year"] %in% d$year, ncol(MD)+1 ] = d$density
  colnames(MD)[ncol(MD)] = substr( tf, 1, nchar(tf)-4)
  
  
} # for all treeFiles


######################################
#' computes the distances for each putative start year of the sample within the chronology
#' @param chronoMD the max density chronology
#' @param sampleMD the sample's max densities
#' @return the distance for each start year
getDistances <- function ( chronoMD, sampleMD ) {
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
# get ranks in chronology for each sample set
######################################

for (l in c(5,10,15,20)) {
  
  
  # get all sample sets
  sampleFiles = list.files( pathSamples , pattern=paste("samples-len-",l,"-set-\\d+\\.csv",sep="") )
  
  allranks = c()
  
  # generate chronology for the set if not already present
  for (f in sampleFiles) {
    
    # read data
    chronoFile = paste(substr(f,1,nchar(f)-4),"chronoMaxDens","csv",sep=".")
    chrono = as.data.frame(read.csv( paste(pathSamples,chronoFile,sep="/"), header=T ))
    # read samples
    samples = as.data.frame(read.csv( paste(pathSamples,f,sep="/") ))
    # compute rank per sample
    sRank = rep(NA, nrow(samples))
    for ( s in 1:nrow(samples) ) {
      # get data of sample
      sMD = MD[ MD$year %in% samples$start[s]:samples$end[s], colnames(MD) %in% samples$file[s] ]
      # compute dist for each shift in chronology
      cDist = getDistances( chrono$maxDens, sMD )
      # get rank for each shift
      cRank = rank(cDist)
      # store rank of sample
      sRank[s] = cRank[ chrono$year %in% samples$start[s]  ]
    } # for each sample
    
    # store ranks of samples
    allranks = cbind(allranks,sRank)
    
  } # for each sampleFile
 
  
  #sum(allranks==1)/sum(!is.na(allranks))
  print(paste("&",l,
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
  #summary(allranks)
  
} # for sample length
