
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

#biomarker Kriterium 2
biomarker.kriterium2 = log2FC.treated.untreated
biomarker.kriterium2[biomarker.kriterium2<0]=-1
biomarker.kriterium2[biomarker.kriterium2>0]=1
biomarker.kriterum2.sum = sapply(1:13299 , function(k) {sum(biomarker.kriterium2[k,]==-1)})
biomarker.kriterium2.sum.names =cbind(rownames(treated), biomarker.kriterum2.sum)
biomarker.kriterium2.sum.names.klein = biomarker.kriterium2.sum.names[,2]<205  #25%
biomarker.kriterium2.sum.names.gross = biomarker.kriterium2.sum.names[,2]>614
komplett.kriterium2 = cbind(biomarker.kriterium2.sum.names , biomarker.kriterium2.sum.names.klein , biomarker.kriterium2.sum.names.gross)
komplett.kriterium2 = komplett.kriterium2[-which(komplett.kriterium2[,4]=="FALSE" & komplett.kriterium2[,3]=="FALSE"),]
