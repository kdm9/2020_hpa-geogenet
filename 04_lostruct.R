#' ---
#' title: Lostruct
#' date: 2021-02-24
#' author: Kevin Murray
#' ---
#'

if (!require("xfun"))              { install.packages("xfun")                                   ; require("xfun")               }
if (!require("tidyverse"))         { install.packages("tidyverse")                              ; require("tidyverse")          }
if (!require("foreach"))           { install.packages("foreach")                                ; require("foreach")            }
if (!require("doParallel"))        { install.packages("doParallel")                             ; require("doParallel")         }
if (!require("parallel"))          { install.packages("parallel")                               ; require("parallel")           }
if (!require("ggplot2"))           { install.packages("ggplot2")                                ; require("ggplot2")            }
if (!require("ggrepel"))           { install.packages("ggrepel")                                ; require("ggrepel")            }
if (!require("ggmap"))             { install.packages("ggmap")                                  ; require("ggmap")              }
if (!require("RColorBrewer"))      { install.packages("RColorBrewer")                           ; require("RColorBrewer")       }
if (!require("SNPRelate"))         { BiocManager::install("SNPRelate")                          ; require("SNPRelate")          }
if (!require("lostruct"))          { remotes::install_github("petrelharp/local_pca/lostruct")   ; require("lostruct")           }
if (!require("windowlickr"))       { remotes::install_github("kdm9/windowlickr")                ; require("windowlickr")        }

knitr::opts_chunk$set(
  warnings=F
)
if (file.exists("data/cache/01_hpa-geogenetics.Rda")) load(file="data/cache/01_hpa-geogenetics.Rda")
theme_set(theme_bw())

NCPUS = as.integer(Sys.getenv("NCPUS", parallel::detectCores(logical=F)))
registerDoParallel(cores=NCPUS)

#' ## BCF windower
#' 
#' So at time of writing, lostruct's vcf parsing code is borked. So, I've
#' written a custom bit of code using windowlickr to do the snp extraction on
#' windows.

windowlickr_windowfun <- function (file, windows=NULL, size=20000, samples=NULL, ...) {
    if (is.null(samples)) { samples = windowlickr:::bcf_getSamples(file) }
    if (is.null(windows)) { windows = windowlickr:::bcf_getWindows(file, windowsize=size, slide=size)}
    pos.fn <- function(n, ...) {
        windows[n,] %>%
            dplyr::select(chrom=contig, start, end=stop)
    }
    win.fn <- function (n,...) {
        if (n > nrow(windows)) stop(paste0("No such window: ", n))
        region = windows[n, "region"]
        ret = windowlickr:::readBCFQuery_(file, region, samples)
        nsnp = length(ret$POS)
        if (nsnp < 1) {
          return(NULL)
        }
        GT = matrix(unlist(ret$GT, recursive = F), nrow=nsnp, byrow = T)
        return(GT)
    }
    attr(win.fn,"max.n") <- nrow(windows)
    attr(win.fn,"region") <- pos.fn
    attr(win.fn,"samples") <- samples
    class(win.fn) <- c("winfun", "function")
    return(win.fn)
}

#' ## Run eigen_windows
#'
#' This is the local PCAs per genome window.

# params: window size & input file
winsize=200000
file="data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.bcf"
npc=20

# Get genome windows and windower function for eigen_windows
windows = windowlickr:::bcf_getWindows(file, windowsize=winsize, slide=winsize)
wf = windowlickr_windowfun(file=file, windows=windows, size=winsize, samples=samp.within.eur)

# Calculate pop. struct. (eigenvectors) of each window
eigwin = xfun::cache_rds({
    eigen_windows(wf, k=npc, mc.cores=NCPUS)
}, file="04_lostruct_eigwin.rds", dir="data/cache/")

# calculate inter-window distances

pcd = xfun::cache_rds({
    pc_dist(eigwin, npc=npc)
}, file="04_lostruct_pcdist.rds", dir="data/cache/")

# exclude NA windows
na.window.idx = is.na(eigwin[,1])
non.na.windows = windows[!na.window.idx,]

# mds
mds =  xfun::cache_rds({
    cmdscale( pcd[!na.window.idx,!na.window.idx], k=5, eig=T)
}, file="04_lostruct_cmdscale.rds", dir="data/cache/")

# make data frame
pts = mds$points
colnames(pts) = paste0("MDS",seq_len(ncol(pts)))
mds.dat = bind_cols(non.na.windows, as.data.frame(pts))
str(mds.dat)

ctgs = windowlickr:::bcf_getContigs(file)
ctg.lens=cumsum(c(0, ctgs$lengths[-length(ctgs$lengths)]))
names(ctg.lens) = ctgs$names

mds.plt = mds.dat %>%
    gather("mds.dim", "value", starts_with("MDS")) %>%
    mutate(overall_pos = ctg.lens[contig] + start)


colours = rep(c("#1f78b4", "#a6cee3"), length.out=nrow(ctgs))
shapes = rep(c(16, 15), length.out=nrow(ctgs))
ggplot(mds.plt, aes(x=overall_pos, y=value)) +
    geom_point(aes(colour=contig, shape=contig)) +
    facet_grid(mds.dim~., space="free_x") +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=shapes) +
    labs(x="Genome position", y="MDS axis") +
    theme(legend.position="none", panel.grid=element_blank())
ggsave("out/plot/04_lostruct-manhattan.svg", width=10, height=6, units="in")
ggsave("out/plot/04_lostruct-manhattan.png", width=10, height=6, units="in", dpi=600)

asai = read_tsv("data/annotation/asaiEffectorBed.txt")

asai.match = asai %>% 
  mutate(contig_matched = sprintf("%06dF|arrow", contig - 1)) %>%
  mutate(overall_pos = ctg.lens[contig_matched] + start)

asai.atr = asai.match %>%
  filter(grepl("ATR", asaiEffector, ignore.case = T))

ggplot(mds.plt, aes(x=overall_pos, y=value)) +
    geom_point(aes(colour=contig, shape=contig)) +
    geom_rug(aes(x=overall_pos, y=0), sides="b", data=asai.match) +
    facet_grid(mds.dim~., space="free_x") +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=shapes) +
    labs(x="Genome position", y="MDS axis") +
    theme(legend.position="none", panel.grid=element_blank())
ggsave("out/plot/04_lostruct-manhattan-asai-rug.png", width=10, height=6, units="in", dpi=600)

ggplot(mds.plt, aes(x=overall_pos, y=value)) +
    geom_point(aes(colour=contig, shape=contig)) +
    geom_rug(aes(x=overall_pos, y=0), sides="b", data=asai.atr) +
    facet_grid(mds.dim~., space="free_x") +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=shapes) +
    labs(x="Genome position", y="MDS axis") +
    theme(legend.position="none", panel.grid=element_blank())
ggsave("out/plot/04_lostruct-manhattan-asai-atr-rug.png", width=10, height=6, units="in", dpi=600)


