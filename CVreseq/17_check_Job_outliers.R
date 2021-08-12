name1 = "Masked_NoInbred_FST_Outflank.igv"
DT1 = read.delim(name1, header = TRUE, sep='\t')
length(DT1$Masked_NoInbred_FST_Outflank)
length(DT1$Masked_NoInbred_FST_Outflank[DT1$Masked_NoInbred_FST_Outflank == 10])

name2 = "Masked_Wild_FST_Outflank.igv"
DT2 = read.delim(name2, header = TRUE, sep='\t')
length(DT2$Masked_Wild_FST_Outflank)
length(DT2$Masked_Wild_FST_Outflank[DT2$Masked_Wild_FST_Outflank == 10])
