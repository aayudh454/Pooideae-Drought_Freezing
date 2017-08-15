res <- res[order(res$padj),]
head(res)
dim(res)
x.sub <- subset(res, padj < 0.01)
head(x.sub)
sub_data=as.data.frame(x.sub)
write.csv(sub_data, file = "Brachypodium_ControlvsDrought.csv", row.names = T, quote = F)
