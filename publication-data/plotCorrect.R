

d = as.data.frame(matrix( 
dimnames= list( NULL,
c("length","ring width","max. density","profile consensus","profile set","two step (MD+PS)")),
c(
  5, 6.6, 15.8, 19.4, 32.2, 42.8
  ,
  10, 28, 44.8, 33.6, 58, 69.6
  ,
  15, 55.4, 74.2, 57.6, 74, NA
  ,
  20, 67.6, 84.8, NA, 82.4, NA
  ,
  25, 75.4, 86.5, NA, NA, NA
),ncol=6, byrow=TRUE))
distCols = !( colnames(d) %in% "length")

matplot(d$length, d[,distCols], lwd=2, type="b", ylim=c(0,100),pch=as.character(1:sum(distCols)),
        ylab="% correct", xlab="sample length (consecutive years)")

legend("bottomright", legend=colnames(d)[distCols], col = 1:5,
       pch=as.character(1:sum(distCols)),
       inset=0.02,
       title="distance:", title.adj = 0, bty="o", bg="white", box.col="white")
