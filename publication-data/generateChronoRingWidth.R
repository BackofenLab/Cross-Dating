
rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathMicaMeta = "../input/PCAB_Ostalb_GD-ring_widths"

# path where sample data is stored
pathSamples = "./samples"

########################################################


######################################
# read and store ring width data
######################################

yearMin = 1916
yearMax = 2004
RW = data.frame(yearMin:yearMax)
colnames(RW) = c("year")

# get all tree files
treeFiles = list.files( pathMicaMeta, pattern=".*\\.meta")
# get available ring width data from each file
for( tf in treeFiles ) {

	d = read.csv(paste(pathMicaMeta,tf,sep="/"))

	RW[ RW[,"year"] %in% d$year, ncol(RW)+1 ] = d$meanNumPoints
	colnames(RW)[ncol(RW)] = substr( tf, 1, nchar(tf)-5)


} # for all treeFiles


######################################
# generate chronology for each sample set
######################################

# get all sample sets
sampleFiles = list.files( pathSamples , pattern="samples-len-\\d+-set-\\d+\\.csv" )

# generate chronology for the set if not already present
for (f in sampleFiles ) {
	# generate name of file to create
	chronoFile = paste(substr(f,1,nchar(f)-4),"chronoRingWidth","csv",sep=".")
	if (file.exists(paste(pathSamples,chronoFile,sep="/"))){ next } # skip if already existing

	# read samples
	samples = read.csv(paste(pathSamples,f,sep="/"))

	# copy and prune data for chronology generation
	cRW = data.frame(RW)
	for (s in 1:nrow(samples)) {
		cRW[ cRW[,"year"] %in% samples[s,"start"]:samples[s,"end"], colnames(RW) %in% as.character(samples[s,"file"]) ] = NA
	} # for samples

	# compute chronology information
	chrono = cbind( as.matrix(cRW[,"year"]), apply(as.matrix(cRW[,2:ncol(cRW)]), 1, mean, na.rm=TRUE) )
	colnames(chrono) = c("year","ringWidth")

	write.csv(chrono, paste(pathSamples,chronoFile,sep="/"),
			quote=FALSE,row.names=FALSE)

} # for all sample sets

######################################
# inspect chronologies
######################################

op <- par(ask=TRUE)
for ( l in c(5,10,15,20) ) {
	d = c();
	for ( f in list.files( pathSamples , pattern=paste("samples-len-",l,"-set-\\d+\\.chronoRingWidth.csv",sep="") ) ) {
		d = cbind( d, read.csv(paste(pathSamples,f,sep="/"))$ringWidth)
	}
	matplot( yearMin:yearMax, d, type="l", main=paste("sample length",l), xlab="year", ylab="ringWidth")
}
par(op)



