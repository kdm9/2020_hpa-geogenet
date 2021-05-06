# ---
# title: Exporatory genetics of HpA
# author: K.D. Murray
# ---

if (!require("tidyverse"))    { install.packages("tidyverse")     ; require("tidyverse")    }
if (!require("foreach"))      { install.packages("foreach")       ; require("foreach")      }
if (!require("doParallel"))   { install.packages("doParallel")    ; require("doParallel")   }
if (!require("parallel"))     { install.packages("parallel")      ; require("parallel")     }
if (!require("ggplot2"))      { install.packages("ggplot2")       ; require("ggplot2")      }
if (!require("ggmap"))        { install.packages("ggmap")         ; require("ggmap")        }
if (!require("fossil"))       { install.packages("fossil")        ; require("fossil")       }
if (!require("vegan"))        { install.packages("vegan")         ; require("vegan")        }
if (!require("gdm"))          { install.packages("gdm")           ; require("gdm")          }
if (!require("ecodist"))      { install.packages("ecodist")       ; require("ecodist")      }
if (!require("SNPRelate"))    { BiocManager::install("SNPRelate") ; require("SNPRelate")    }
if (!dir.exists("out/plot")) dir.create("out/plot", recursive=T)


# # Genotype data
#
# We will use SNPRelate for most genetics. First, we need to create a file in
# SNPrelate's GDS format. This takes ages so we cache the result, and only
# create a new gds if it's missing.

if (!file.exists("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds")) {
    cat("Updating GDS file from VCF\n")
    snpgdsVCF2GDS("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.vcf.gz",
                  "data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds")
} else {
    cat("Not updating gds file\n")
}

gds = snpgdsOpen("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds", allow.duplicate=T)
gds.sum  = snpgdsSummary(gds)

# ## Broad data exploration: missingness, freqs
#
# Before doing anything interesting with the data, we should inspect the SNP
# matrix. Specifically, we wish to thin out samples that have too much missing
# data, and also remove hyper-rare SNPs and those with large amounts of missing
# data.

# ### Sample missingness

samp.miss = snpgdsSampMissRate(gds)
hist(samp.miss, main="Sample missing prop")

# Looks broadly good, about three crappy samples with < 80% completeness. Will go
# ahead with those in there though as they still look reasonable. 
# Setting the threshold at 20% seems reasonable.

sample.missing.thresh = 0.2
abline(v=sample.missing.thresh)


# ### SNP missingness
#
# To get both the MAF/AF/missingness per SNPs we use snpgdsSNPRateFreq. 

snp.rate = snpgdsSNPRateFreq(gds)

# Let's see what the missingness rate looks like across SNPs.

hist(snp.rate$MissingRate, main="SNP missing prop")
snp.missing.thresh = 0.3
abline(v=snp.missing.thresh)

# Per-snp missingness is very good. Will set the snp missing threshold at 30%
# as that removes the worst few percent of SNPs

mean(snp.rate$MissingRate>0.30)

# ### SNP MAF

hist(snp.rate$AlleleFreq, main="SNP AF", breaks=50)
hist(snp.rate$MinorFreq, main="SNP MAF", breaks=50)

# Fairly nice folded SFS, but it seems the reference allele is frequently the
# rarer allele. There aren't many hyper-rare alleles (typically there is a
# massive spike around 1%), so we won't filter on MAF for basic pca/distance
# based analyses.
#
# ### Create combined eu/geno sample list
#
# In the previous notebook we got the list of Eurasian samples. We want to
# combine this with out geneitc filtering here to pass both those subsettings
# forward to future analyses.

samp.geno.ok = gds.sum$sample.id[samp.miss < sample.missing.thresh]
writeLines(samp.geno.ok, "data/metadata/samples_genotype_ok.txt")
samp.within.eur = readLines("data/metadata/samples_within_europe.txt")
samp.within.eur.geno.ok = intersect(samp.within.eur, samp.geno.ok)

# # Distance-based Popgen
#
# Now for the very first *genetic* analysis: distance-based exploration of
# population structure. Specifically, we will use both direct PCA-based
# methods, and heirarchical clustering of IBS distances.

# ## IBS dendrogram
#
# Identity by state is a simple genetic distance: the distance between two
# samples is the hamming distance between their chromosomes. We first run this
# on all samples, including genetic failures, then on just those samples within
# Eurasia that pass the missing data filter above. In each case we use only
# those SNPs that passed the missing data threshold we set above (which is ~97%
# of the data).

snp.ibs.all = snpgdsIBS(gds, missing.rate = snp.missing.thresh, num.thread = 4)

snp.ibs.all %>%
    snpgdsHCluster() %>%
    snpgdsCutTree() %>%
    snpgdsDrawTree()

# So in this set of all samples there appear to be either two or three major
# groups, depending on where you draw the line. There are a few samples that
# appear weird, and a reasonable number of close relatives as expected in this
# largely population-based sampling design.

snp.ibs.goodeu = snpgdsIBS(gds, missing.rate = snp.missing.thresh,
                           sample.id=samp.within.eur.geno.ok, num.thread = 4)

snp.ibs.goodeu %>%
    snpgdsHCluster() %>%
    snpgdsCutTree() %>%
    snpgdsDrawTree()

## When removing the american and poor quality samples, not much changes in the
## dendrogram. We can also plot the distance matrix directly:

image(1- snp.ibs.goodeu$ibs)

# ## PCA

snp.pca = snpgdsPCA(gds, missing.rate = snp.missing.thresh)
plot(snp.pca)

# ## Align to (EU) geno data
#
# Now match the rows/cols of the IBS matrix to get our distance

m = match(samp.within.eur,  snp.ibs$sample.id)
geno.dist = 1 - snp.ibs$ibs[m, m]
dim(geno.dist)
length(geo$ind)
rownames(geno.dist) = colnames(geno.dist) = samp.within.eur
hist(geno.dist)

# ## Simple IBD

geo.dist = as.matrix(eu.geo[, c("longitude", "latitude")]) %>%
    earth.dist() %>%
    as.dist()
ibd.geno.dist = geno.dist %>%
    as.dist()
vegan::mantel(ibd.geno.dist, geo.dist)
mg = ecodist::mgram(ibd.geno.dist, geo.dist)
mg$mgram[,4]
plot(mg, xlab="Distance (km)")


# ## GDM


gdm.geno = cbind(ind=rownames(geno.dist), as.data.frame(as.matrix(geno.dist)))
gdm.preds = eu.geo %>%
    dplyr::select(ind, latitude, longitude)
gdm.geodist = cbind(ind=eu.geo$ind, as.data.frame(as.matrix(geo.dist)))
dim(gdm.geno)
dim(gdm.preds)
dim(gdm.geodist)


sp = formatsitepair(gdm.geno, 3, predData=gdm.preds, siteColumn="ind",
                    XColumn="longitude", YColumn="latitude",
                    distPreds=list("geography.gcd"=gdm.geodist)) %>%
    filter(s2.matrix_1 > 5) # only include site pairs with a geographic distance > 5km
mdl = gdm(sp, geo=F)
plot(mdl, include.rug=T, rug.sitepair=sp)
