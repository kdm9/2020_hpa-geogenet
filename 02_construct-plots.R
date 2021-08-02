install.packages(c("ggspatial", "ggplot2", "sf", "rnaturalearth",
                   "rnaturalearthdata", "ggforce"))

library(tidyverse)
library(ggspatial)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggforce)
library(foreach)
theme_set(theme_bw())

meta = read_tsv("data/metadata/europe-metadata-geoclust.tsv")
load("data/cache/02_construct_admix_prop.Rds", verbose=T)

# # Pie map plots

plt  = crossval_admix_prop %>%
    select(mdl, K, rep, both) %>%
    unnest(both) %>%
    mutate(group = dist(cbind(longitude, latitude)) %>%
                       hclust() %>%
                       cutree(h=0.8)) %>%
    group_by(mdl, K, rep, group, pop) %>%
    summarise(latitude=mean(latitude), longitude=mean(longitude),
              proportion=mean(proportion), size=n()) %>%
    ungroup() %>%
    group_by(mdl, K, rep, group) %>%
    nest() %>%
    mutate(plot=purrr::map(data,
               function(d)
                    geom_arc_bar(aes(x0=longitude, y0=latitude, r0=0,
                                     fill=pop,r=0.5, amount=proportion),
                                data=d, stat="pie", inherit.aes=F, linetype="blank")))


eu.bbox = c(left=-12, right=29, top=60, bottom=34)
eu = ne_countries(scale = "medium", returnclass = "sf", continent="europe")

ggplot(eu) +
    geom_sf(fill="white", colour="grey") +
    plt %>%
       filter(mdl=="sp", K==2, rep==1) %>%
       pull(plot) +
    #scale_colour_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") +
    coord_sf(xlim = eu.bbox[1:2], ylim = eu.bbox[c(4,3)], expand = FALSE) +
    labs(x=NULL, y=NULL) +
    theme(legend.position="none", panel.grid=element_blank())
ggsave("out/Hpa_cs_v2/Hpa_cs_k2-piemap.png", height=6.3, width=7, unit="in", dpi=1200)

ggplot(eu) +
    geom_sf(fill="white", colour="grey") +
    plt %>%
       filter(mdl=="sp", K==3, rep==1) %>%
       pull(plot) +
    #scale_colour_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") +
    coord_sf(xlim = eu.bbox[1:2], ylim = eu.bbox[c(4,3)], expand = FALSE) +
    labs(x=NULL, y=NULL) +
    theme(legend.position="none", panel.grid=element_blank())
ggsave("out/Hpa_cs_v2/Hpa_cs_k3-piemap.png", height=6.3, width=7, unit="in", dpi=1200)

# # Structure bar plots

k23_cs = crossval_admix_prop %>%
    filter(K %in% 2:3, rep==1) 


samp.within.eur.geno.ok = readLines("data/metadata/samples_within_europe_genotype_ok.txt")
admix.indiv = readLines("data/admixture/individualOrdered.txt")

k23_admix = foreach(K=2:3, rep=1:10, .combine=rbind) %do%
{
    admix = read.table(sprintf("data/admixture/rep%d_P/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.vcf.gz.%d.Q", rep, K), header=F) %>%
        mutate(indiv=admix.indiv) %>%
        tidyr::gather("pop", "admix.prop", -indiv) %>%
        mutate(pop = as.integer(sub("^V", "" ,pop)), rep=rep, K=K)

}


admix.meta = meta %>%
    filter(ind %in% samp.within.eur.geno.ok) %>%
    select(-pop)

k23_admix_pop = k23_admix %>% 
    inner_join(admix.meta, by=c("indiv"="ind")) %>%
    group_by(pop, K, rep, geo.clust) %>%
    summarise_at(vars(admix.prop, latitude, longitude), mean) %>%
    mutate(set=sprintf("ADMIXTURE K%d", K))

k23_cs_pre = k23_cs %>%
    select(-size, -longitude, -latitude) %>%
    mutate(set=sprintf("conStruct K%d %s", K, c("sp"="Spatial", "nsp"="Non-spatial")[mdl]), pop=as.integer(pop)) %>%
    rename(admix.prop=proportion)


barplot.dat = bind_rows(k23_cs_pre, k23_admix_pop)

kview(k23_cs)
kview(barplot.dat)