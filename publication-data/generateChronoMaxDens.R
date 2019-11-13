
rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathMaxDens = "../input/PCAB_Ostalb_GD-max_densities"

# path where sample data is stored
pathSamples = "./samples"

########################################################


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
# generate chronology for each sample set
######################################

# get all sample sets
sampleFiles = list.files( pathSamples , pattern="samples-len-\\d+-set-\\d+\\.csv" )

# generate chronology for the set if not already present
for (f in sampleFiles ) {
	# generate name of file to create
	chronoFile = paste(substr(f,1,nchar(f)-4),"chronoMaxDens","csv",sep=".")
	if (file.exists(chronoFile)){ next } # skip if already existing

	# read samples
	samples = read.csv(paste(pathSamples,f,sep="/"))

	# copy and prune data for chronology generation
	cMD = data.frame(MD)
	for (s in 1:nrow(samples)) {
		cMD[ cMD[,"year"] %in% samples[s,"start"]:samples[s,"end"], colnames(MD) %in% as.character(samples[s,"file"]) ] = NA
	} # for samples

	# compute chronology information
	chrono = cbind( as.matrix(cMD[,"year"]), apply(as.matrix(cMD[,2:ncol(cMD)]), 1, mean, na.rm=TRUE) )
	colnames(chrono) = c("year","maxDens")

	write.csv(chrono, paste(pathSamples,chronoFile,sep="/"),
			quote=FALSE,row.names=FALSE)

} # for all sample sets

######################################
# inspect chronologies
######################################

op <- par(ask=TRUE)
for ( l in c(5,10,15,20) ) {
	d = c();
	for ( f in list.files( pathSamples , pattern=paste("samples-len-",l,"-set-\\d+\\.chronoMaxDens.csv",sep="") ) ) {
		d = cbind( d, read.csv(paste(pathSamples,f,sep="/"))$maxDens )
	}
	matplot( yearMin:yearMax, d, type="l", main=paste("sample length",l), xlab="year", ylab="maxDens")
}
par(op)



