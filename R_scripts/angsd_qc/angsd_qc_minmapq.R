mafs = c()
for (x in 10:37) {
  file = paste0("test_", x, ".mafs")
  df <- read.delim(file, header = TRUE, sep='\t')
  #print(df$knownEM[1])
  maf <- df$knownEM[2]
  mafs <- append(mafs, maf)
}

plot(x=seq(10, 37, 1), y=mafs)
mafs[20]-mafs[10]


