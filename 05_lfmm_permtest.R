# # LFMM permutation tests

library(tidyverse)
library(LEA)
library(SNPRelate)
library(foreach)
library(doParallel)
if (!dir.exists("data/cache/05_lfmm/")) dir.create("data/cache/05_lfmm")
load("data/cache/05_lfmm/lfmm-inputs.Rds", verbose=T)

Sys.setenv("OMP_NUM_THREADS"=1)
Sys.setenv("OPENBLAS_NUM_THREADS"=1)
Sys.setenv("MKL_NUM_THREADS"=1)

NCPUS = as.integer(Sys.getenv("NCPUS", Sys.getenv("NSLOTS", parallel::detectCores(logical=F))))
registerDoParallel(cores=NCPUS)
cat(paste("Using", NCPUS, "cores\n"))

# ## Permutation tests
#
# Let's form an expected null the dumb way: permute the environment data rows
# (**not** each var independently) to break any statstical link between
# genotypes and the environment.

perms = xfun::cache_rds({
    foreach(i=1:100, .combine=bind_rows) %dopar% {
        cat(sprintf("replicate %d\n", i))
        perm.env = env4lea[sample(seq_len(nrow(env4lea))),]
        cat(sprintf("permenv %d\n", i))

        rep.lf = lfmm2("data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                       perm.env, K=nmf.K)

        cat(sprintf("lfmm2 %d\n", i))
        rep.pv = lfmm2.test(rep.lf, "data/tmp/lfmm/genotypes.lfmm_imputed.lfmm",
                            perm.env, linear = TRUE)
        cat(sprintf("test %d\n", i))

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
}, file="permutationtests", dir="data/cache/05_lfmm/")
