# ---
# title: Construct on Hpa
# date: 2021-01-02
# author: Kevin Murray
# ---

# Construct models population structure as a combination of discrete and
# continuous structure.

if (!require("tidyverse"))         { install.packages("tidyverse")                      ; require("tidyverse")          }
if (!require("foreach"))           { install.packages("foreach")                        ; require("foreach")            }
if (!require("doParallel"))        { install.packages("doParallel")                     ; require("doParallel")         }
if (!require("parallel"))          { install.packages("parallel")                       ; require("parallel")           }
if (!require("ggplot2"))           { install.packages("ggplot2")                        ; require("ggplot2")            }
if (!require("ggrepel"))           { install.packages("ggrepel")                        ; require("ggrepel")            }
if (!require("ggmap"))             { install.packages("ggmap")                          ; require("ggmap")              }
if (!require("RColorBrewer"))      { install.packages("RColorBrewer")                   ; require("RColorBrewer")       }
if (!require("fossil"))            { install.packages("fossil")                         ; require("fossil")             }
if (!require("vegan"))             { install.packages("vegan")                          ; require("vegan")              }
if (!require("SNPRelate"))         { BiocManager::install("SNPRelate")                  ; require("SNPRelate")          }
if (!require("conStruct"))         { install.packages("conStruct")                      ; require("conStruct")          }
if (!require("constructhelpers"))  { remotes::install_github("kdm9/constructhelpers")   ; require("constructhelpers")   }


if (!dir.exists("out/Hpa_cs_v2/")) dir.create("out/Hpa_cs_v2/")

NCPUS = as.integer(Sys.getenv("NCPUS", parallel::detectCores(logical=F)))
registerDoParallel(cores=NCPUS)
cat(paste("Using", NCPUS, "cores\n"))


# ## Geographic clustering
#
# First, we need to select samples and find geographic clusters.

meta = read_tsv("data/metadata/europe-metadata.tsv")
geo.dist.indiv = meta %>%
    column_to_rownames(var = "ind") %>%
    dplyr::select(longitude, latitude) %>%
    fossil::earth.dist()
geo.clust.indiv = geo.dist.indiv %>%
    hclust() %>%
    cutree(h=5) # cut at 5km radius
geo.cluster.name = sprintf("GC%02d", geo.clust.indiv)
meta$geo.clust = geo.cluster.name
write_tsv(meta, "data/metadata/europe-metadata-geoclust.tsv")

samp.within.eur.geno.ok = readLines("data/metadata/samples_within_europe_genotype_ok.txt")

# So at a population radius of 5km, we have `r length(geo.cluster.name)`
# clusters.
#
# Here are the cluster sizes (y clusters of size x): 

table(table(geo.cluster.name))

# ## Basic genetics
#
# Here we need to extract the SNP data, and reduce it to allele frequencies at
# each geographic cluster, and remove snps with null variance.

# First,  get the indivs by snps matrix (122x500k)

gds = snpgdsOpen("data/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.gds", allow.duplicate=T)
gds.sum  = snpgdsSummary(gds)
gen = snpgdsGetGeno(gds, sample.id = samp.within.eur.geno.ok, with.id = T)

# Now, we reduce the 122x500k matrix of SNPs to a 68x500k matrix of population
# allele frequencies. We skip columns (SNPs) with null variance. 

gen.pop = xfun::cache_rds({
    N.SNP = ncol(gen$genotype)
    foreach(i=seq_len(N.SNP), .combine=cbind) %do% {
        v = var(gen$genotype[, i], na.rm=T)
        if (!is.finite(v) || v <=0) {
            return(NULL)
        }
        tapply(gen$genotype[, i], geo.cluster.name, function(x) {mean(x, na.rm=T)/2})
    }
}, file="02_01_genpop", dir="data/cache/",  compress="xz")


# ## Run construct
#
# Here we do some data prep, then run all the construct runs needed to do cross
# validation (NB: we don't use construct's own cross validation here as it runs
# all K sequentially which makes it take K times as long as it needs to.


geo.clust.dat = meta %>%
    group_by(geo.clust, pop) %>%
    summarise(latitude = mean(latitude),
              longitude = mean(longitude)) %>%
    ungroup()
geo.clust.latlong = geo.clust.dat %>%
    dplyr::select(longitude, latitude) %>% 
    as.matrix()
geo.clust.dist = geo.clust.latlong %>%
    fossil::earth.dist() %>%
    as.matrix()


# So to run construct we shell out to the cluster, as it takes ages. So, we save
# the inputs to an Rds and then load the results again here for plotting and
# interpretation.

cxv = xfun::cache_rds({
    constructhelpers::csh.x.validation(
        prefix="out/Hpa_cs_v2/HpA_cs_v2",
        freqs=gen.pop,
        coords=geo.clust.latlong,
        geoDist=geo.clust.dist,
        K=1:8,
        train.prop=0.9,
        n.reps=16,
        n.iter=20000,
        n.nodes=NCPUS,
        save.files=F,
        make.figs=T)
}, file="02_02_cs_xval", dir="data/cache/",  compress="xz")


# # Summarise construct cross validation

# Here we plot the log likelihoods for each run. NB that in theory, the 

cxv.std.fit = cxv %>%
    select(-construct.res, -mdl.fit) %>%
    group_by(mdl) %>%
    mutate(std.fit = mean.fit - max(mean.fit)) %>%
    ungroup()
ggplot(cxv.std.fit, aes(K, std.fit)) +
    geom_point(position="jitter") +
    facet_wrap(~mdl)+
    theme_bw()
ggsave("out/Hpa_cs_v2/Hpa_cs_xval_ll.pdf", width=8,height=5)


layer = cxv %>%
    group_by(mdl, K) %>%
    mutate(layer.contrib =  purrr::map2(construct.res, data.block, function (x, y) {
            lc = conStruct:::calculate.layer.contribution(x,y)
            lr = 1:length(lc)
            tibble(layer = lr, layer.contrib=sort(unlist(lc), dec=T))
        })
    ) %>%
    select(mdl, K, rep, layer.contrib) %>%
    unnest(layer.contrib)

plot.dat =  layer %>%
    filter(K>1) %>%
    mutate(
        K=as.factor(sprintf("K=%d", K)),
        layer=as.factor(layer),
        rep=as.factor(sprintf("Rep%02d", rep)),
        mdl=fct_relevel(ifelse(mdl == "sp", "Spatial", "Non-Spatial"), "Spatial", "Non-Spatial")
    )
ggplot(plot.dat, aes(x=layer, y=layer.contrib)) +
    geom_violin() +
    #geom_point(aes(colour=rep)) + 
    #scale_colour_manual(values=rainbow(length(levels(plot.dat$rep))), name="Replicate") +
    labs(y="Layer Contribution", x="Layer (ancestral popn.)") +
    facet_grid(mdl ~ K, scales="free_x", space = "free_x") +
    theme_bw()
ggsave("out/Hpa_cs_v2/Hpa_cs_xval_lc.pdf", width=8,height=5)


# ## K=2 and K=3 Piecharts
#
# Here we use conStruct's `make.all.the.plots` to generate plots for each run for K in 1:4.


if (!dir.exists("out/Hpa_cs_v2/matp/")) dir.create("out/Hpa_cs_v2/matp/", recursive=T)

cxv %>%
    filter(K %in% 1:4, rep %in% 1:4) %>%
    mutate(code=sprintf("K%d_%s_rep%d", K, mdl, rep)) %>%
    select(code, construct.res, data.block) %>%
    purrr::pmap(function(code, construct.res, data.block) {
            conStruct::make.all.the.plots(list(chain_1=construct.res), data.block, paste0("out/Hpa_cs_v2/matp/Hpa_cs_v2_", code, ".pdf"))
        }) %>%
    invisible()


# # Paper plots
#
# We need two plots for the paper: a structure-style barplot with both Gautam's
# ADMIXTURE result, and the barplots from conStruct, and a pie-chart figure but
# overlaid on a basemap of europe.

crossval_admix_prop = cxv %>%
    filter(mdl == "sp", K %in% 2:3, rep==1) %>%
    mutate(coords = purrr::map(data.block,
               function(x) 
                   x[["coords"]]%>%
                       as_tibble() %>%
                       mutate(site=1:nrow(.))),
           admix.prop = purrr::map(construct.res,
               function(x)
                   x[["MAP"]][["admix.proportions"]] %>%
                       magrittr::set_colnames(sprintf("pop_%d", 1:ncol(.))) %>%
                       as_tibble()),
           both=purrr::map2(coords, admix.prop,
               function(x, y)
                   bind_cols(x, y) %>%
                       pivot_longer(cols=starts_with("pop_"),
                                    names_to="pop",
                                    names_prefix="pop_",
                                    values_to="proportion"))
           ) %>%
    select(mdl, K, rep, both) %>%
    unnest(both) %>%
    mutate(group = dist(cbind(longitude, latitude)) %>%
                       hclust() %>%
                       cutree(h=0.8)) %>%
    group_by(mdl, K, rep, group, pop) %>%
    summarise(latitude=mean(latitude), longitude=mean(longitude),
              proportion=mean(proportion), size=n()) %>%
    ungroup()

save(crossval_admix_prop, file="data/cache/02_construct_admix_prop.Rds")
