#' ---
#' title: Halfmaximal LD over Hpa
#' date: 2021-01-06
#' author: Kevin Murray
#' ---

if (!require("boringLD"))   { remotes::install_github("kdm9/boringLD") ; require("boringLD")   }
if (!require("tidyverse"))  { install.packages("tidyverse")            ; require("tidyverse")  }
if (!require("foreach"))    { install.packages("foreach")              ; require("foreach")    }
if (!require("doParallel")) { install.packages("doParallel")           ; require("doParallel") }
if (!require("parallel"))   { install.packages("parallel")             ; require("parallel")   }
if (!require("ggplot2"))    { install.packages("ggplot2")              ; require("ggplot2")    }

NCPUS = as.integer(Sys.getenv("NCPUS", parallel::detectCores(logical=F)))
registerDoParallel(cores=NCPUS)
cat(paste("Using", NCPUS, "cores\n"))
theme_set(theme_bw())
if (!dir.exists("out")) dir.create("out")
if (!dir.exists("data/cache")) dir.create("data/cache")
if (!dir.exists("out/LD")) dir.create("out/LD")

samp.within.eur = readLines("data/metadata/samples_within_europe_genotype_ok.txt")

#' # Calculate LD across genome

#' Here we calculate the distance to half-maximal decay of $R^2$
#' in 100kbp windows (tiled every 50kb) across the HpA genome.
#' This is a big computation that we cache to make life easier.


halfmax = xfun::cache_rds({
    windowed_halfmax(bcf="data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.bcf",
                     windowsize=100000,
                     minMAF=0.1,
                     maxMissing=0.4)
}, file="03_01_hpa-halfmax", dir="data/cache/",  compress="xz")

#' Now plot this over the genome.

halfmax.plot = halfmax %>%
    arrange(contig, start) %>%
    mutate(culmpos = cumsum(stop - start))

ggplot(halfmax.plot, aes(culmpos, halfmax)) +
    geom_line() +
    labs(x="Genome Position", y="Half-maximal R2 distance (bp)") +
    scale_y_log10()
ggsave("out/LD/halfmax-LD-over-genome.pdf", width=8, height=5)

#' So, by eye, LD extends between 10kb and 100kb, varying
#' considerably across the genome.

#' ## Pairwise corr plot for a contig
#'
#' Here we want a heatmap style plot of correlations among SNPs
#' on one contig. First, we need to get the genotypes

ctg = windowlickr:::bcf_getContigs("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.bcf")
geno = windowlickr:::bcf_getGT(
    "data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.bcf",
    region=sprintf("%s:1-%d", ctg$names[1], ctg$lengths[1]),
    minMAF=0.1,
    maxMissing=0.4,
    rowsAreSamples=T,
    samples=samp.within.eur
)

#' As making all pairwise r2s takes too much ram, we limit it to
#' 10k snps. The first of these blocks of code is randomly
#' sampled over the first contig.

transparent.black = rgb(0, 0, 0, alpha = 15, maxColorValue = 255)
which.snp.randsubset = sort(sample(seq_len(geno$nSNP))[1:10000])
minigt = geno$GT[,which.snp.randsubset]
correl = cor(minigt, use="pairwise.complete.obs") ^ 2
bpdist = dist(geno$POS[which.snp.randsubset])
image(correl, main="Subsampled ctg1")
plot(correl ~ as.matrix(bpdist), pch=".", col=transparent.black)

#' This next block is the first 10k snps (about 1/3 the contig)
#' without subsampling.

which.snp.ld = 1:10000
minigt = geno$GT[,which.snp.ld]
correl = cor(minigt, use="pairwise.complete.obs") ^ 2
bpdist = dist(geno$POS[which.snp.ld])
image(correl, main="first10k ctg1")
plot(correl ~ as.matrix(bpdist), pch=".", col=transparent.black)
pdf("out/LD/first10ksnp.pdf")
image(correl, main="first10k ctg1")
plot(correl ~ as.matrix(bpdist), pch=".", col=transparent.black)
dev.off()
