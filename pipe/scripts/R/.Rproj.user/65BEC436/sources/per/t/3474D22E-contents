# Can do This for Maize - Affy -------------------------------------------------
# Read mappings, values per probe
vals <- read.table(gzfile("~/tmp/CDF/ZM/ZM-00014/TXT/no_controls/GSM684664.rma.gz"),
                   header=T, stringsAsFactors = F)
map <- read.table("~/tmp/CDF/ZM/ZM-00014/from72/best_genelmeasure_ZM_AFFY_MAIZE_Gene.csv", 
                  header=F, stringsAsFactors = F, skip = 15, sep = ",")

# map between the two
vals$probeset_norm <- tolower(gsub("_at$", "", vals$probeset))
map[, 1] <- tolower(map[, 1])
vals$probeset_map <- map$V2[match(vals$probeset_norm, map$V1)]

# drop nonexisting maps and duplicates
vals_test <- vals[!is.na(vals$probeset_map), c("probeset_map", "signal", "stderror")]
vals_test <- vals_test[which(!(duplicated(vals_test$probeset_map) | rev(duplicated(vals_test$probeset_map)))), ]
# read "truth", values attributed to probesets by old pipeline (old cdf)
vals_truth <- read.table(gzfile("~/tmp/CDF/ZM/ZM-00014/from72/GSM684664.rma.gz"), header=T, stringsAsFactors = F)

# susbet to common values, order one by the other
rownames(vals_test) <- vals_test$probeset_map
rownames(vals_truth) <- vals_truth$probeset
combi_vec <- match(rownames(vals_test), rownames(vals_truth))
vals_truth_compare <- vals_truth[combi_vec[!is.na(combi_vec)], ]
assertthat::are_equal(sum(rownames(vals_truth_compare) == rownames(vals_test)), nrow(vals_test))

signal_diff <- vals_truth_compare$signal - vals_test$signal
summary(signal_diff)
summary(vals_test$signal)
summary(vals_truth_compare$signal)
plot(signal_diff)

################################################################################
# Can do this for Agilent as well ----------------------------------------------
vals <- read.table(gzfile("~/tmp/CDF/MM/MM-00534/TXT/GSM1360685_US90503634_252800514459_S01_GE2_107_Sep09_1_1.rma.gz"),
                   header=T, stringsAsFactors = F)
map <- read.table("~/tmp/CDF/MM/MM-00534/from72/best_genelmeasure_MM_AGIL_8x60K_Gene.csv", 
                  header=F, stringsAsFactors = F, skip = 15, sep = ",")

# map between the two
vals$probeset_norm <- tolower(gsub("_at$", "", vals$probeset))
map[, 1] <- tolower(map[, 1])
vals$probeset_map <- map$V2[match(vals$probeset_norm, map$V1)]

# drop nonexisting maps and duplicates
vals_test <- vals[!is.na(vals$probeset_map), c("probeset_map", "signal", "stderror")]
vals_test <- vals_test[which(!(duplicated(vals_test$probeset_map) | rev(duplicated(vals_test$probeset_map)))), ]
# read "truth", values attributed to probesets by old pipeline (old cdf)
vals_truth <- read.table(gzfile("~/tmp/CDF/MM/MM-00534/from72/GSM1360685_US90503634_252800514459_S01_GE2_107_Sep09_1_1.rma.gz"), header=T, stringsAsFactors = F)

# susbet to common values, order one by the other
rownames(vals_test) <- vals_test$probeset_map
rownames(vals_truth) <- vals_truth$probeset
combi_vec <- match(rownames(vals_test), rownames(vals_truth))
vals_truth_compare <- vals_truth[combi_vec[!is.na(combi_vec)], ]
assertthat::are_equal(sum(rownames(vals_truth_compare) == rownames(vals_test)), nrow(vals_test))

signal_diff <- vals_truth_compare$signal - vals_test$signal
summary(signal_diff)
summary(vals_test$signal)
summary(vals_truth_compare$signal)
plot(signal_diff)

################################################################################
# Summary Affy: ----------------------------------------------------------------

# > summary(signal_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5202.2   254.6   720.6  1436.2  1731.3 36554.4 
# > summary(vals_test$signal)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 35.22   433.48  1114.82  2129.65  2583.92 23985.59 
# > summary(vals_truth_compare$signal)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 52.32   740.93  1891.94  3565.89  4315.08 42359.24

# lost only 11009 - 10441 = 568 genes


# Summary Agilent: -------------------------------------------------------------

# > summary(signal_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -841695       0       0    4707      25 3970309 
# > summary(vals_test$signal)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27     388     876   22129    6066 3445306 
# > summary(vals_truth_compare$signal)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27     403    1116   26837    8868 3972947 

# lost only 22334 - 21326 = 1008 genes

################################################################################
# Comparing DE from GV ---------------------------------------------------------

old <- read.table("~/tmp/CDF/MM/MM-00534/fromGV/DE_AgilentCDF_FDR01.csv",
                  header=T, stringsAsFactors = F, sep = ",")
# old <- subset(old, select = -c(Description, p.value))
new <- read.table("~/tmp/CDF/MM/MM-00534/fromGV/DE_NebionCDF_FDR01.csv",
                  header=T, stringsAsFactors = F, sep = ",")
# new <- subset(new, select = -c(Description, p.value))

rownames(old) <- old$Gene
rownames(new) <- new$Gene

# susbet to common values, order one by the other
idxer <- rownames(new) %in% rownames(old)
new_yes <- new[idxer, ]
new_no <- new[!idxer, ] # 27 Genes missing from 569, good
new_no[, c("Log.ratio", "minFDR", "Gene.Symbol")]

idxer <- rownames(old) %in% rownames(new)
old_yes <- old[idxer, ]
# old_new <- new[!idxer, ]
old_yes <- old_yes[match(rownames(new_yes), rownames(old_yes)), ]

assertthat::are_equal(sum(rownames(old_yes) == rownames(new_yes)), nrow(new_yes))

plot(old_yes$minFDR, type = "p", col = "blue", xlab="index", ylab="minFDR")
title("FDR for commonnly appearing samples")
points(new_yes$minFDR, col = "orange", xlab="", ylab="")
legend("topright", legend = c("agilent", "custom"), col = c("blue", "orange"), 
       pch = c("o", "o"))

plot(old_yes$Log.ratio, type = "p", col = "blue", xlab="index", ylab="Log.ratio")
title("Log.ratio for commonnly appearing samples")
points(new_yes$Log.ratio, col = "orange", xlab="", ylab="")
legend("topright", legend = c("agilent", "custom"), col = c("blue", "orange"), 
       pch = c("o", "o"))

