# # Construct-related plots for the paper

#install.packages(c("ggspatial", "ggplot2", "sf", "rnaturalearth",
#                   "rnaturalearthdata", "ggforce", "rgeos"))
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
load("data/cache/02_construct_admix_prop-aaz.Rds", verbose=T)

csxvap = crossval_admix_prop %>%
    select(mdl, K, rep, both) %>%
    unnest(both) %>%
    mutate(geo.clust = sprintf("GC%02d", site))


# # Pie map plots

plt  = csxvap %>%
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


eu.bbox = c(left=-12, right=55, top=60, bottom=34)
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
ggsave("out/plot/02_Hpa_cs-aaz_k2-piemap.svg", height=6.3, width=7, unit="in", dpi=1200)

ggplot() +
    plt %>%
       filter(mdl=="sp", K==2, rep==1) %>%
       pull(plot) +
    coord_fixed()
ggsave("out/plot/02_Hpa_cs-aaz_k2-justpies.svg", height=6.3, width=7, unit="in", dpi=1200)

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
ggsave("out/plot/02_Hpa_cs-aaz_k3-piemap.svg", height=6.3, width=7, unit="in", dpi=1200)

ggplot() +
    plt %>%
       filter(mdl=="sp", K==3, rep==1) %>%
       pull(plot) +
    coord_fixed()
ggsave("out/plot/02_Hpa_cs-aaz_k3-justpies.svg", height=6.3, width=7, unit="in", dpi=1200)

# # Structure bar plots

samp.within.eur.geno.ok = readLines("data/metadata/samples_within_europe_genotype_ok.txt")

k23_cs = csxvap %>%
    filter(K %in% 2:3, rep==1) 


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
    mutate(set=sprintf("AM K%d", K))

k23_cs_pre = k23_cs %>%
    mutate(set=sprintf("CS K%d %s", K, mdl), pop=as.integer(pop)) %>%
    rename(admix.prop=proportion)


barplot.dat = bind_rows(k23_cs_pre, k23_admix_pop) %>%
    select(set, K, rep, geo.clust, pop, admix.prop) %>%
    mutate(pop=as.factor(pop))

ggplot(barplot.dat, aes(x=geo.clust, y=admix.prop)) +
    geom_bar(aes(colour=pop, fill=pop), stat="identity", position="stack") +
    facet_grid(.~set) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_classic() +
    labs(x=NULL, y=NULL) +
    theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(),
        strip.text.x=element_text(angle=90))
ggsave("out/plot/02_construct-admixture-barplots.svg", width=3, height=6)
