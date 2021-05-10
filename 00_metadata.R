# ---
# title: Basic metadata processing of HpA data
# author: K.D. Murray
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

library(tidyverse)
library(raster)
library(ggmap)
library(ggplot2)
library(GGally)
theme_set(theme_bw())
if (!dir.exists("data/metadata")) dir.create("data/metadata", recursive=T)
if (!dir.exists("out/plot")) dir.create("out/plot", recursive=T)


# # Geographic exploration
#
# ## Read metadata
#
# This metadata is on all samples, including the few american ones. This is
# direct from Gautam, and includes the fix for the issue with the swapped
# coordinates for the Irish samples (2-3 samples from Ireland had lat/longs in
# Spain due to an ID number confusion).

geo = read_csv("data/finalInds.AllGeoInfo.csv")
summary(geo[, c("latitude", "longitude")])


# ## Maps
#
# This plots samples for just the european samples on a map.

eu.bbox = c(left=-18, right=32, top=62, bottom=32)
basemap.terr = get_stamenmap(eu.bbox, zoom=5, maptype="terrain-background")
eu.geo = geo %>%
    filter(latitude > eu.bbox["bottom"], latitude < eu.bbox["top"],
           longitude > eu.bbox["left"], longitude < eu.bbox["right"])


ggmap(basemap.terr, darken=c(0.4, "white"), legend="topleft", extent="device") +
    geom_point(aes(longitude, latitude, colour=pop), data=eu.geo,
               alpha=1, size=2) +
    labs(x="Longitude", y="Latitude") +
    scale_colour_discrete(name="Population") + 
    theme_bw()
ggsave("out/plot/00_samples-within-europe-map.svg")

# To pass this forwards to all future analyses we save a list of samples within
# Eurasia. We are focused on the Eurasian samples here as this is the native
# range of both Hpa and Ath

samp.within.eur = eu.geo %>%
    pull(ind)
write_tsv(eu.geo, "data/metadata/europe-metadata.tsv")
writeLines(samp.within.eur, "data/metadata/samples_within_europe.txt")


# # Climate
#
# In order to associate our G with an E in a GEA, we need some 'E'.
# Specifcically, Bioclim variables for each population/sampling point.
# We need to use the raster package to extract these from the WorldClim2
# Bioclim 30s layer TIFFs. 

if (!file.exists("data/metadata/europe-metadata-env.tsv")) {
    rasters = do.call("stack", lapply(1:19,
        function (i) raster(sprintf("/data/kevin/work/gis/WorldClim2/wc2.0_bio_30s_%02d.tif", i),
                varname=sprintf("Bio%02d", i))))

    meta.env = eu.geo %>%
        dplyr::filter(ind %in% samp.within.eur)
    coordinates(meta.env) = ~ longitude + latitude

    meta.env =  bind_cols(
        meta.env %>%
        as.data.frame(),
        extract(rasters, meta.env) %>%
        as.data.frame() %>%
        rename_with(function(x) gsub("wc2.0_bio_30s_([0-9]+)", "Bio\\1", x, perl=T))
    )

    write_tsv(meta.env, "data/metadata/europe-metadata-env.tsv")
} else {
    meta.env = read_tsv("data/metadata/europe-metadata-env.tsv")
}

# Let's plot the 19 bioclim variables and lat/long to see any overall trends.
# We expect some significant inter-variable correlation, as many of the bioclim
# variables are essentially reformulations of each other.

meta.env %>%
    dplyr::select(-ind, -loc, -pop) %>%
    as.data.frame() %>%
    ggpairs(upper=list(continuous="points", combo="facethist", discrete="facetbar", na="na"))
ggsave("out/plot/00_env-variable-pairs.png", height=20, width=20)
