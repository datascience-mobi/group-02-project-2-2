
#Specific Analysis

#Find Biomarker for cisplatin through FC
cisplatin.col=c()
j=1

while(j<820)
{
  if(isTRUE(meta[j,3]== "cisplatin"))
  {cisplatin.col= c(cisplatin.col,j)
  }
  j = j +1
}

FC.cisplatin = log2FC.treated.untreated[,cisplatin.col]
mean.FC.rows <- apply(FC.cisplatin, 1, mean)
sorted.mean = sort(mean.FC, decreasing = TRUE)

#finding/visualizing most extreme FC values for cisplatin
lowest.FC <- sorted.mean[13290:13299]
par(mar = c(5, 5, 5, 5))
barplot(lowest.FC, horiz = TRUE, xlim = c(-1.3,0), main= "lowest log2 FC-values for cisplatin", xlab= "mean log2FC values in different celllines", col= "firebrick", names.arg=c("FTO", "PLK1", "VPS8", "POLR38", "MAPKAP1", "STAG1", "TBCD", "C11orf49","ANKS1A", "COMMD10"), las=1, cex.names =0.8)

highest.FC <- sorted.mean[1:10]
par(mar = c(5, 5, 5, 5))
barplot(highest.FC, horiz = TRUE, xlim = c(0,1.8), main= "highest log2 FC-values for cisplatin", xlab= "mean log2FC values in different celllines", col= "lightgreen",las=1, cex.names =0.8)

#upordownregulation
upordown = (treated.scaled/untreated.scaled)
upordown = upordown[,c(95:149)]
upordown = t(upordown)
upordown.biomarker = upordown[ , c("SUPT3H","BBS9", "ANKRA2", "GNBL1", "KDM4C", "DDIT3", "TTC28","NBEA", "PCCA","PPP1R15A", "FTO", "PLK1", "VPS8", "POLR3B", "MAPKAP1", "STAG1", "TBCD", "C11orf49", "ANKS1A", "COMMD10")]
upordown.biomarker = t(upordown.biomarker)
heatmap(upordown.biomarker, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

#Kriterium 2
is.neg = FC.cisplatin<0
i =1
j=1
a=1
biomarker2.row = c()
while(j<13300){
  while(i<56){
    if(is.neg[j,i]==TRUE)
     {a=a+1}
    i=i+1}
  if(a>49)
     {biomarker2.row= c(biomarker2.row, j)}
  if(a<6)
    {biomarker2.row= c(biomarker2.row, j)}
  a=1
  i=1
  j=j+1}
#das sind jetzt die nrow von Genen dich sich unter cisplatin
#in 50 von 55 Fällen in die gleiche Richtung verändern
biomarker2 = c(rownames(FC.cisplatin[biomarker2.row]))