dir="/lab/data/seqcore/5124-NM/10x_analysis_5124-NM"

summaryStats = list()
for (i in 1:8) {
  filename = file.path(dir, paste0("Sample_5124-NM-", i), "summary.csv")
  print(filename)
  summaryStats[[i]] = read.csv(filename)
}

ss = do.call("rbind", summaryStats)
write.table(ss, "/lab/work/ccrober/T1D_U01/5124-NM.summary.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
