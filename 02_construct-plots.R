
install.packages(c("ggspatial", "ggplot2", "sf", "rnaturalearth",
                   "rnaturalearthdata", "ggforce"))
library(ggspatial)
library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggforce)
theme_set(theme_bw())

plt  = crossval_admix_prop %>%
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
