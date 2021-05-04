# # LFMM GEA
#
# This notebook runs an LFMM2-based Genotype-environment Association over the
# SNPs derived from Gautam's variant calls.

library(tidyverse)
library(LEA)
library(raster)
library(SNPRelate)

load("data/sample_selection.Rda", verbose=T)


# # Metadata

meta = read_csv("data/finalInds.AllGeoInfo.csv")


# # Climate
#
# In order to associate our G with an E in a GEA, we need some 'E'.
# Specifcically, Bioclim variables for each population/sampling point.
# We need to use the raster package to extract these from the WorldClim2
# Bioclim 30s layer TIFFs. 

rasters = do.call("stack", lapply(1:19,
    function (i) raster(sprintf("/data/kevin/work/gis/WorldClim2/wc2.0_bio_30s_%02d.tif", i),
                        varname=sprintf("Bio%02d", i))))

meta.env = meta %>%
    dplyr::filter(ind %in% samp.within.eur)
coordinates(meta.env) = ~ longitude + latitude

meta.env =  bind_cols(
    meta.env %>%
        as.data.frame(),
    extract(rasters, meta.env) %>%
        as.data.frame() %>%
        rename_with(function(x) gsub("wc2.0_bio_30s_([0-9]+)", "Bio\\1", x, perl=T))
)

# Let's plot the 19 bioclim variables and lat/long to see any overall trends.
# We expect some significant inter-variable correlation, as many of the bioclim
# variables are essentially reformulations of each other.

meta.env %>%
    dplyr::select(-ind, -loc, -pop) %>%
    as.data.frame() %>%
    plot()


# # Genetics
#
# The G for the GEA. We will use Gautam's variant calls

gds = snpgdsOpen("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds")

snp.stats = snpgdsSNPRateFreq(gds, sample.id=samp.within.eur)

# ## SNP Filtering
#
# We want a very high quality set of SNPs, so our SNP filtering thresholds will
# be generally more severe than for distance-based analyses like PCA or GDMs.
#
# Whatever their functional relevance, SNPs with low MAF are highly unlikely to
# be statstically significantly associated with the environment with only
# 100-ish indiviuals. We filter out all SNPs with < 10% MAF, i.e. present in 13
# samples, as this does not reduce the number of SNPs too much, and with 13
# samples we have some statsitcal power. Note that for a "real" GEA hit, we
# would expect nearby, common variation to "tag" causal variation as it does in
# a GWAS. The extent to which this occurs is governed by the number of
# alleles/haplotypes involved in each indiviual tradjectory of adaption to the
# environment. 

snp.maf.thresh = 0.1
hist(snp.stats$MinorFreq)
abline(v=snp.maf.thresh)

# Similarly, we filter quite severely on missingness too. This is particularly
# important as LFMM2 requires a complete genotype matrix, and we will used NNMF
# to impute the missing genotypes. Luckily we have excellent matrix
# comletelness already, and lose only a small number of SNPs in the long tail
# of sparser data.

snp.missing.thresh = 0.1
hist(snp.stats$MissingRate)
abline(v=snp.missing.thresh)


snp.mask = snp.stats$MinorFreq >= snp.maf.thresh & snp.stats$MissingRate < snp.missing.thresh
mean(snp.mask)
snp.mask.id = which(snp.mask)


gn.list = snpgdsGetGeno(gds, sample.id=samp.within.eur, snp.id=snp.mask.id, with.id=T, verbose=T)
gn = gn.list$genotype
rownames(gn) = gn.list$sample.id
colnames(gn) = gn.list$snp.id
dim(gn)

# Althogether, we keep about 74% of SNPs, which is plenty. We can now move on
# to imputing missing genotypes.

# ## Imputation
#
# LFMM requires a complete genotype data matrix. Luckily the LEA package
# provides a basic form of genotype imputation based on non-negative matrix
# factorisation (NNMF). We will apply imputation, and assess the correlation of
# imputed genotypes to ground truth.


# Unfortunately, the LEA package is a bit bonkers in terms of data flow, in
# that one must write out each matrix or part thereof to independent files. It
# even requires that we write the data to different formats depending on the
# function used. Therefore, we will save the genotype data to both "geno" and
# "lfmm" formats, and the environment variables to the "env" format.


# first, set the missing values to 9 as required
gn4lea = gn
gn4lea[is.na(gn4lea)] = 9

# then write genotype data to both a .geno and a .lfmm file
if (!dir.exists("data/lfmm")) dir.create("data/lfmm")
write.geno(gn4lea, "data/lfmm/genotypes.geno")
write.lfmm(gn4lea, "data/lfmm/genotypes.lfmm")

# and then env data
env4lea = meta.env %>%
    dplyr::select(-loc, -pop) %>%
    remove_rownames() %>%
    column_to_rownames("ind")
write.env(env4lea, "data/lfmm/env-latlon.env")


# ### PCA
#
# A precursor to imputation is PCA

pc = pca("data/lfmm/genotypes.lfmm", scale=T)
tw = tracy.widom(pc)
plot(10 * -log10(tw$pvalues))
which(tw$pvalues < 1e-4)


nmf = snmf("data/lfmm/genotypes.geno",
           K=1:30,
           CPU=14,
           iterations=100,
           entropy=T,
           repetitions=10)

plot(nmf)
