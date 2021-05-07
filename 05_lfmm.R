# # LFMM GEA
#
# This notebook runs an LFMM2-based Genotype-environment Association over the
# SNPs derived from Gautam's variant calls.

library(tidyverse)
library(LEA)
library(SNPRelate)
library(foreach)
library(doParallel)
if (!dir.exists("data/cache/05_lfmm/")) dir.create("data/cache/05_lfmm")
load("data/sample_selection.Rda", verbose=T)


# # Metadata

meta = read_csv("data/finalInds.AllGeoInfo.csv")
samp.within.eur.geno.ok = readLines("data/metadata/samples_within_europe_genotype_ok.txt")

# In order to associate our G with an E in a GEA, we need some 'E'.
# Specifcically, Bioclim variables for each population/sampling point.
# We need to use the raster package to extract these from the WorldClim2
# Bioclim 30s layer TIFFs. 
#
# I have added the climate variables to the metadata in a previous notebook,
# the results of which I load below.

meta.env = read_tsv("data/metadata/europe-metadata-env.tsv")

# # Genetics
#
# The G for the GEA. We will use Gautam's variant calls

gds = snpgdsOpen("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds")

snp.stats = snpgdsSNPRateFreq(gds, sample.id=samp.within.eur.geno.ok)

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
# to imputing missing genotypes. Finally, get some per-snp metadata like
# positions etc.

contigs = windowlickr:::bcf_getContigs("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.vcf.gz") %>%
    mutate(offsets = cumsum(c(0, lengths[-length(lengths)])))
snp.list = snpgdsSNPList(gds, sample.id=samp.within.eur)[snp.mask,]


# # LFMM
#
# We wish to run an LFMM2 (Caye et al. 2019). There are many intermediate steps.

# ## Data Munging


# Unfortunately, the LEA package is a bit bonkers in terms of data flow, in
# that one must write out each matrix or part thereof to independent files. It
# even requires that we write the data to different formats depending on the
# function used. Therefore, we will save the genotype data to both "geno" and
# "lfmm" formats, and the environment variables to the "env" format.


# first, set the missing values to 9 as required
gn4lea = gn
gn4lea[is.na(gn4lea)] = 9

# then write genotype data to both a .geno and a .lfmm file
if (!dir.exists("data/tmp/lfmm")) dir.create("data/tmp/lfmm")
write.geno(gn4lea, "data/tmp/lfmm/genotypes.geno")
write.lfmm(gn4lea, "data/tmp/lfmm/genotypes.lfmm")

# and then env data
env4lea = meta.env %>%
    dplyr::select(-loc, -pop) %>%
    remove_rownames() %>%
    column_to_rownames("ind")
write.env(env4lea, "data/tmp/lfmm/env-latlon.env")


# ## PCA
#
# A precursor to imputation is PCA

pc = pca("data/tmp/lfmm/genotypes.lfmm", scale=T)
tw = tracy.widom(pc)
plot(10 * -log10(tw$pvalues))
which(tw$pvalues < 1e-4)


# ## SNMF
#
# Non-negative matrix factorisation decomposes a matrix (our genotypes) into
# the product of two matrices (P and Q). P and Q describe respectively the allele
# frequencies of a series of ancestral sub-populations, and the admixture
# proportions of each current indiviual with resepect to these ancestral
# populations. This is in essence the same problem that STRUCTURE and ADMIXTURE
# address.


nmf = xfun::cache_rds({
    snmf("data/tmp/lfmm/genotypes.geno",
               K=1:20,
               CPU=14,
               iterations=100,
               entropy=T,
               repetitions=4)
}, file="snmf.Rds", dir="data/cache/05_lfmm/")

# Now we plot the cross-entropy (the SNMF loss function).  We are running SNMF
# on our population as the imputation is done based on this model of admixture.
# Therefore, unlike when using construct or ADMIXTURE, we care more about the
# statsitcal properties of the data than the biological reality of any detected
# ancestral populations. In other words, in order to more accurately impute
# genotypes, we want this model of ancestral struture to err on the side of
# overfitting, using a higher value of K than might seem biogicaly realistic if
# it is statstically supported.

plot(nmf)

# The above plot of cross-entropy reaches its minimum at approximately K=18, so
# we will use that as our value of K. We also need to pick the best run across
# our replicates.

nmf.K = 18
nmf.best = which.min(cross.entropy(nmf, K = nmf.K))


# We may as well also plot a barchart, but as the value of K is unrealistically
# high for any biological inference, it won't make a whole bunch of sense.

barchart(nmf, K = nmf.K, run = nmf.best)


# ## Imputation
#
# LFMM requires a complete genotype data matrix. Luckily the LEA package
# provides a basic form of genotype imputation based on non-negative matrix
# factorisation (NNMF). We will apply imputation, and assess the correlation of
# imputed genotypes to ground truth.

impute(nmf, "data/tmp/lfmm/genotypes.lfmm",
       method="mode", K=nmf.K, run=nmf.best)

# TODO: drop out 1% of snps, assess accuracy


# ## LFMM
#
# Now we run the model at the core of LFMM, a mixed model that associates G ~ E
# + LF. This has two parts, first the fitting of the latent factors, and then a
# simple linear scan test (`lfmm2.test`).

lfmm.res = lfmm2("data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                 "data/tmp/lfmm/env-latlon.env",
                 K=nmf.K)

lfmm.pv = lfmm2.test(lfmm.res, "data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                     "data/tmp/lfmm/env-latlon.env",  linear = TRUE)

# In order to plot a manhattan plot, we need some per-snp metadata, as well the
# pvalues. We need to calculate both the fdr, and for plotting purposes -log10(p).

snp.dat = snp.list %>%
    left_join(contigs, by=c("chromosome"="names")) %>%
    dplyr::transmute(snp.id, chromosome, position, overall.pos=position + offsets)

plot.dat = lfmm.pv$pvalue %>%
    t() %>%
    as.data.frame() %>%
    magrittr::set_colnames(colnames(env4lea)) %>%
    bind_cols(snp.dat, .) %>%
    pivot_longer(names_to="envvar", values_to="pval", -c(1:4)) %>%
    mutate(
        fdr=p.adjust(pval, method="BH"),
        log.fdr=-log10(fdr)
    )

# And now the plots. First, one massive manhattan plot.

ggplot(plot.dat, aes(x=overall.pos, y=log.fdr)) +
    geom_point(aes(colour=chromosome)) +
    scale_colour_manual(values=rep(c("#a6cee3", "#b2df8a"), length=nrow(contigs)), guide=F) +
    theme_bw() +
    facet_grid(envvar~.) +
    labs(y="-log10(fdr)", x="Genome Position (bp)")
ggsave("out/plot/05_lfmm-manhattan-plots.svg", width=300, height=210, units="mm")
ggsave("out/plot/05_lfmm-manhattan-plots.png", width=300, height=210, units="mm")


# Then plot a histogram of pvalues for each bioclim variable.

ggplot(plot.dat, aes(x=pval)) +
    geom_histogram() +
    theme_bw() +
    facet_wrap(~envvar, ncol=3)
ggsave("out/plot/05_lfmm-pval-hist.svg", width=20, height=30, units="cm")
ggsave("out/plot/05_lfmm-pval-hist.png", width=20, height=30, units="cm")


# And some Q-Q plots of P-values. There is some pretty full on pvalue inflation happening

qq.dat = plot.dat %>%
    dplyr::select(envvar, pval) %>%
    group_by(envvar) %>%
    arrange(pval) %>%
    mutate(expected=ppoints(pval))

ggplot(qq.dat, aes(x=-log10(expected), y=-log10(pval))) +
    geom_point() +
    theme_bw() +
    facet_wrap(~envvar, ncol=3, scales="free")
#ggsave("out/plot/05_lfmm-qqplot.svg", width=20, height=30, units="cm")
ggsave("out/plot/05_lfmm-qqplot.png", width=20, height=30, units="cm")


# ## Permutation tests
#
# Let's form an expected null the dumb way: permute the environment data rows
# (**not** each var independently) to break any statstical link between
# genotypes and the environment.

perms = xfun::cache_rds({
    foreach(i=1:100, .combine=bind_rows) %do% {
        cat(sprintf("replicate %d\n", i))
        perm.env = env4lea[sample(seq_len(nrow(env4lea))),]

        rep.lf = lfmm2("data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                       perm.env, K=nmf.K)
        rep.pv = lfmm2.test(rep.lf, "data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                            perm.env, linear = TRUE)

        rep.pv$pvalue %>%
            t() %>%
            as.data.frame() %>%
            magrittr::set_colnames(colnames(perm.env)) %>%
            bind_cols(snp.dat, .) %>%
            pivot_longer(names_to="envvar", values_to="pval", -c(1:4)) %>%
            mutate(
                fdr=p.adjust(pval, method="BH"),
                log.fdr=-log10(fdr),
                rep=i,
            )
    }
}, file="permutationtests.Rda", dir="data/cache/05_lfmm/")
