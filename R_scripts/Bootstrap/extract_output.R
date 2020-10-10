
ps_005 = c()
ps_001 = c()

for(prefix in seq(100, 2700, 200)){
  filename = paste0(prefix, '_grab.size.csv')
  print(filename)
  dat = read.delim(filename, header = TRUE, sep=',')
  ps_005 = c(ps_005, sum(dat$adj.ps < 0.05))
  print(ps_005)
  ps_001 = c(ps_001, sum(dat$adj.ps < 0.01))
}

plot(seq(100, 2700,200), ps_005)
plot(seq(100, 2700,200), ps_001)
