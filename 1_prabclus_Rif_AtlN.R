library(prabclus)

## test Rif vs N-Atl

setwd("~/Work/Acantho/IBD/praclust_method/1_Rif_vs_N-Atl/")

# Import geographic coordinates and calculate distances

geographic_coordinates <- read.csv("Coordinates_AtlN_Rif.csv", sep="\t", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./Rif_NAtl_genet_dist.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./Rif_NAtl_genet_dist.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
#genetic_distances <- as.dist(genetic_distances, diag = F, upper = F)
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples
#write.table(genetic_distances, "./Genet_dist_Rif_NAtl.csv", sep=";")

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="9_loci_Alleles_N_Atl_v_Rif.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_AtlN_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_AtlN <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_AtlN_a ~ geo_dist_within_AtlN, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_a ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_AtlN_a) ~ c(geo_dist_within_AtlN)))
abline(lm(c(genet_dist_within_Rif_a) ~ c(geo_dist_within_Rif)), col="red")
points(genet_dist_between ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between) ~ c(geo_dist_between)), col="blue")

# Create plots and do reg tests using shared alleles distances

genet_dist_within_AtlN_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_AtlN <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- H01_reg$pval[2]*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H031_reg$pval, "\n", "p-value H032 ", H032_reg$pval, sep="")}

plot(genet_dist_within_AtlN_s ~ geo_dist_within_AtlN, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_s ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_AtlN_s) ~ c(geo_dist_within_AtlN)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_Rif_s) ~ c(geo_dist_within_Rif)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)
