
rm( list = ls() )

### IN/OUT SETUP #######################################

# path where MICA alignment meta information is stored
pathMicaMeta = "../input/PCAB_Ostalb_GD-ring_widths"

# path where to store output (overwriting existing files)
pathOutput = "./samples"

sampleTarget = matrix(c(
10, 50, 5
,
15, 33, 7
,
20, 25, 10
), byrow=TRUE, ncol=3)
colnames(sampleTarget) = c("len","num","numSets")

########################################################

# get all tree files
treeFiles = list.files( pathMicaMeta, pattern=".*\\.meta")

# create output path if not existing
dir.create(pathOutput, showWarnings=FALSE)

# for each sample set
for (s in 1:nrow(sampleTarget)) {

###########################################
# aggregate possible sample ranges
###########################################

# container to store possible sample ranges
pool <- matrix(NA, ncol=3,nrow=0)
colnames(pool) <- c("file","start","end")
pool <- as.data.frame(pool)

# current sample length
sLen = sampleTarget[s,"len"];

# get possible samples for each file
for( tf in treeFiles ) {
	# read data
	d <- cbind(read.csv(paste(pathMicaMeta,tf,sep="/"))$year,1)
	# stop if not enough years
	if (nrow(d) < sLen) { next; }
	# identify consecutive years
	for (r in 2:nrow(d)) {
		if (d[r,1] == d[r-1,1]+1) { d[r,2] = d[r-1,2]+1 }
	}
	# store putative sample ends
	x <- as.numeric(d[ d[,2] >= sLen ,1])
	tfPrefix = substr(tf, 1, nchar(tf)-5)
	for ( y in x ) {
		pool[nrow(pool)+1,] = c(tfPrefix, y-sLen+1, y)
	}
} # for treeFiles

for ( setId in 1:sampleTarget[s,"numSets"] ) {

###########################################
# pick random set of non-overlapping samples
###########################################

maxTries = 9999
tryCount = 0
repeat {
tryCount = tryCount +1
if (tryCount > maxTries) { stop(paste("Tried",tryCount,"times to find non-overlapping sample set.. giving up now!")) }
# get random set of samples
sNum = sampleTarget[s,"num"]
samples = sort(sample(1:nrow(pool),sNum,replace=FALSE))
# check if non-overlapping
samplesNonOverlapping = TRUE
if (sNum > 1) {
for ( sidx in 2:length(samples) ) {
	si = samples[sidx]
	sp = samples[sidx-1]
	samplesNonOverlapping = pool[sp,"file"] != pool[si,"file"] || as.numeric(pool[sp,"end"]) < as.numeric(pool[si,"start"])
	if ( ! samplesNonOverlapping ) { break }
}
}
# stop sampling if non-overlapping set was found
if (samplesNonOverlapping) { break }
} # repeat until samples are non-overlapping


###########################################
# store samples
###########################################

write.csv(pool[samples,],
	paste(pathOutput,"/samples-len-",sLen,"-set-",setId,".csv",sep=""), 
	quote=FALSE,row.names=FALSE)

} # for each sample set to generate per length

} # for samples (rows in sampleTarget)



