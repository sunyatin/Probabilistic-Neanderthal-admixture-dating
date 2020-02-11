# Rscript function.R infile minD maxD fileChromWeights

args = commandArgs(trailingOnly=TRUE)
file = args[1]
minD = as.numeric(args[2])
maxD = as.numeric(args[3])
fileChrom = args[4]

print(file)
print(minD)
print(maxD)
print(fileChrom)

weighted_jackknife = function(chrom, values, chrom_weights_matrix) {
  # Busing F, Meijer E, Leeden R (1999) Delete-m Jackknife for Unequal m. Statistics and Computing 9: 3–8.
  
  if (length(chrom)!=length(values)) stop('chrom vector should be the same size as values vector')
  if (ncol(chrom_weights_matrix)!=2) stop('chrom_weights_matrix should be a matrix with two columns: chrom ID, chrom weight')
  
  chrom = as.character(chrom)
  general_mean = mean(values, na.rm=T)
  chrom_means = c(by(values, chrom, FUN=mean, na.rm=T))
  chrom_weights_matrix = data.frame(q = as.character(chrom_weights_matrix[,1]),
                                    w = as.numeric(as.character(chrom_weights_matrix[,2])))
  
  CHR = data.frame(q=names(chrom_means), mean=chrom_means)
  CHR = merge(CHR, chrom_weights_matrix, by='q', sort=F, all.x=TRUE, all.y=FALSE)
  
  mrem = sapply(seq_along(CHR$mean), function(i) mean(CHR$mean[-i], na.rm=T))
  CHR$mean = mrem
  
  n = sum(CHR$w)
  g = nrow(CHR)
  
  mean = g*general_mean - sum( (1-CHR$w/n)*CHR$mean  )
  
  H = n/CHR$w
  term1 = sum( (1-CHR$w/n)*CHR$mean )
  term2 = H*general_mean - (H-1)*CHR$mean - g*general_mean + term1
  var = 1/g * sum( 1/(H-1)*term2**2 )
  sd = sqrt(var)
  
  return(c('general_mean'=general_mean,
           'jackknife_mean'=mean,
           'jackknife_sde'=sd))
}

fit_decay_curve_archaicadmixture = function(X, minD = 0.1, maxD = 20) {
  X = subset(X, x>=minD & x<=maxD)
  
  options(show.error.messages = FALSE)
  coef = c('A'=NA, 't'=NA, 'c'=NA)
  fit = try( nls(y~A*exp(-x*t)+c, data=X, start=list(A=.1, t=1, c=1e-3)), silent = TRUE )
  if (class(fit) != "try-error") coef = coef(fit)
  options(show.error.messages = TRUE)
  
  return(coef)
}


#######
##############
#####################
##############
#######

X = read.table(file, header=T, sep='\t')
CHR = read.table(fileChrom, header=F)

#print('Setting NA covariances as 0')

NAMES = c('## Bin(cM)', 'Covariance')

#ggplot(X) + geom_point(aes(bin.right.bound, cov)) + facet_wrap(~chrom, scales='free')

xy = X
xy$cov = xy$cov * xy$n.pairs
xy = na.omit(xy[,c('bin.right.bound', 'cov', 'n.pairs')])
xy <- aggregate(xy[,c("cov", "n.pairs")], list(xy$bin.right.bound), sum)
xy = xy[order(xy$Group.1),]
xy$y = xy[,2]/xy[,3]
xy <- subset(xy, Group.1>=minD & Group.1<=maxD)
xy = data.frame('#Bin(cM)'=xy$Group.1, 'Covariance'=xy$y)
#xy[is.na(xy$Covariance),2] = 0
write.table(xy, paste0(file,'.hard'), row.names=F, quote=F, col.names = NAMES, sep=' ')

FIT = c()
for (q in 1:22) {
  xx = subset(X, chrom==q & bin.right.bound>=minD & bin.right.bound<=maxD)
  xx = na.omit(xx[,c('bin.right.bound', 'cov')])
  names(xx) = c('x', 'y')
  mod = fit_decay_curve_archaicadmixture(xx, minD, maxD)
  t = mod['t'] * 100
  FIT = rbind(FIT, data.frame(chrom=q, y=t))
}

write.table(FIT, paste0(file,'.hard.jin'), row.names=F, quote=F, sep=' ')

jk = weighted_jackknife(FIT$chrom, FIT$y, CHR)
write.table(jk, paste0(file,'.hard.fit'), row.names=F, quote=F, sep=' ')


################

    
xy = X
xy$weighted.cov = xy$weighted.cov * xy$sum.log10.lhoods
xy = na.omit(xy[,c('bin.right.bound', 'weighted.cov', 'sum.log10.lhoods')])
xy <- aggregate(xy[,c("weighted.cov", "sum.log10.lhoods")], list(xy$bin.right.bound), sum)
xy = xy[order(xy$Group.1),]
xy$y = xy[,2]/xy[,3]
xy <- subset(xy, Group.1>=minD & Group.1<=maxD)
xy = data.frame('#Bin(cM)'=xy$Group.1, 'Covariance'=xy$y)
cat("#", file = paste0(file,'.soft'))
#xy[is.na(xy$Covariance),2] = 0
write.table(xy, paste0(file,'.soft'), row.names=F, quote=F, sep = ' ', col.names = NAMES)

FIT = c()
for (q in 1:22) {
  xx = subset(X, chrom==q & bin.right.bound>=minD & bin.right.bound<=maxD)
  xx = na.omit(xx[,c('bin.right.bound', 'weighted.cov')])
  names(xx) = c('x', 'y')
  mod = fit_decay_curve_archaicadmixture(xx, minD, maxD)
  t = mod['t'] * 100
  FIT = rbind(FIT, data.frame(chrom=q, y=t))
}

write.table(FIT, paste0(file,'.soft.jin'), row.names=F, quote=F, sep=' ')

jk = weighted_jackknife(FIT$chrom, FIT$y, CHR)
write.table(jk, paste0(file,'.soft.fit'), row.names=F, quote=F, sep=' ')



