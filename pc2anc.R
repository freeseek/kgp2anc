#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (c) 2016-2018 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

args <- commandArgs(trailingOnly = TRUE)
pca <- args[1]
kgp <- args[2]
out <- args[3]

# load and merge input tables
df1 <- read.table(pca, header=TRUE)
df2 <- read.table(kgp, header=TRUE)
df <- merge(df1, df2, all=TRUE)

# establish number of PCs to use
if (length(args) < 4) {
  npc <- 5
} else {
  npc <- as.numeric(args[4])
}

df$EAS <- NaN

# European populations
eur <- df$POP %in% c('CEU', 'TSI', 'FIN', 'GBR', 'IBS') | grepl("^EUR", df$POP)
df[which(eur), 'EUR'] <- 1; df[which(eur), 'AFR'] <- 0; df[which(eur), 'NAT'] <- 0; df[which(eur), 'EAS'] <- 0; df[which(eur), 'SAS'] <- 0

# African populations
afr <- df$POP %in% c('YRI', 'LWK', 'GWD', 'MSL', 'ESN') | grepl("^AFR", df$POP) & !grepl("ASW", df$POP) & !grepl("ACB", df$POP)
df[which(afr), 'EUR'] <- 0; df[which(afr), 'AFR'] <- 1; df[which(afr), 'NAT'] <- 0; df[which(afr), 'EAS'] <- 0; df[which(afr), 'SAS'] <- 0

# East Asian populations
eas <- df$POP %in% c('CHB', 'JPT', 'CHS', 'CDX', 'KHV') | grepl("^EAS", df$POP)
df[which(eas), 'EUR'] <- 0; df[which(eas), 'AFR'] <- 0; df[which(eas), 'NAT'] <- 0; df[which(eas), 'EAS'] <- 1; df[which(eas), 'SAS'] <- 0

# South Asian populations
sas <- df$POP %in% c('GIH', 'PJL', 'BEB', 'STU', 'ITU') | grepl("^SAS", df$POP)
df[which(sas), 'EUR'] <- 0; df[which(sas), 'AFR'] <- 0; df[which(sas), 'NAT'] <- 0; df[which(sas), 'EAS'] <- 0; df[which(sas), 'SAS'] <- 1

# Admixed populations
mixed <- (df$POP %in% c('ASW', 'ACB', 'MXL', 'PUR', 'CLM', 'PEL') | grepl("^AMR", df$POP) | grepl("ASW", df$POP) | grepl("ACB", df$POP)) & !(df$IID %in% c('HG01880', 'HG01944'))
for (anc in c('EUR', 'AFR', 'NAT')) {
  df[which(mixed), anc] <- df[which(mixed),anc] * 1/(1-df[which(mixed),"UNK"])
}
df[which(mixed), 'EAS'] <- 0
df[which(mixed), 'SAS'] <- 0

# predict global ancestral proportions
superpops <- c('EUR', 'AFR', 'NAT', 'EAS', 'SAS')
for (anc in superpops) {
  fit <- lm(as.formula(paste0(anc, ' ~ ', paste0('PC', 1:npc, collapse = ' + '))), df)
  idx <- is.na(df[, anc])
  pred <- predict.lm(fit, df[idx,], interval = 'prediction')
  df[idx, anc] <- pred[, 'fit']
  df[idx, paste0(anc, '_L95')] <- pred[, 'lwr']
  df[idx, paste0(anc, '_U95')] <- pred[, 'upr']
}

# output data to disk
idx <- df$POP=='SET' & !is.na(df$POP)
col.names <- c('FID', 'IID', 'SEX', superpops, paste0(superpops, '_L95'), paste0(superpops, '_U95'))
write.table(format(df[idx, col.names], digits = 3, scientific = FALSE), out, sep = '\t', row.names = FALSE, quote = FALSE)
