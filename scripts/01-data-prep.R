# Data preparation script
# Stenio I A Foerster
# 16 Dec 2025
# https://foersterst.github.io


# Traits --------------------------------------------------------------------------------------


# Libraries & data

library(tidyverse)
library(ape)
library(rptR)

# community matrix
comm <- readxl::read_excel("data/scorp-comm.xlsx", sheet = 1, skip = 1)

# species list
sp <- colnames(comm)[grepl("_", colnames(comm))]
sp <- sort(substr(sp, 6, 99))

# raw morphological data
# read in raw data
r1 <- readxl::read_excel("data/scorp-trait-base.xlsx", sheet = 1, skip = 1)
r2 <- readxl::read_excel("data/scorp-trait-base.xlsx", sheet = 2, skip = 1)

# combine both tabs & adjust species names to match tip labels
r1 |>
  dplyr::select(-c(source, specimen)) -> r1
r2 |>
  dplyr::select(colnames(r1)) -> r2

raw_data <- rbind(r1, r2)
rm(r1, r2)

setdiff(sp, raw_data$species)

# Calculate repeatability to know if interspecific trait variation is larger than intraspecific (including sexual dimorphism)
# Resulting R values should be as large as possible.

# complete data
raw_data |>
  filter(species %in% sp) |>
  mutate(
    cheL_rel = cheL / carL,
    cheW_rel = cheW / carL,
    cheD_rel = cheD / carL,
    carL_log10 = log10(carL)
  ) |>
  select(species, sex, carL_log10, cheL_rel, cheW_rel, cheD_rel) |>
  drop_na() -> data_ssd

# repeatability for carapace length
rpt(
  carL_log10 ~ sex + (1 | species),
  grname = "species",
  data = data_ssd,
  datatype = "Gaussian",
  nboot = 100
)

# repeatability for relative chela length
rpt(
  cheL_rel ~ sex + (1 | species),
  grname = "species",
  data = data_ssd,
  datatype = "Gaussian",
  nboot = 100
)

# repeatability for relative chela width
rpt(
  cheW_rel ~ sex + (1 | species),
  grname = "species",
  data = data_ssd,
  datatype = "Gaussian",
  nboot = 100
)

# repeatability for relative chela height
rpt(
  cheD_rel ~ sex + (1 | species),
  grname = "species",
  data = data_ssd,
  datatype = "Gaussian",
  nboot = 100
)

# Conclusion: repeatability is always >50%, meaning that trait variation across species is larger than within species.

# species means
raw_data |>
  filter(species %in% sp) |>
  mutate(
    cheL_rel = cheL / carL,
    cheW_rel = cheW / carL,
    cheD_rel = cheD / carL,
    carL_log10 = log10(carL)
  ) |>
  group_by(species) |>
  summarise(
    carL_log10 = mean(carL_log10, na.rm = T),
    cheL_rel = mean(cheL_rel, na.rm = T),
    cheW_rel = mean(cheW_rel, na.rm = T),
    cheD_rel = mean(cheD_rel, na.rm = T),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA, .x))) -> traits

# read in functional groups
fn <- readxl::read_excel("data/func-groups.xlsx", sheet = 1)

# merge
traits <- full_join(traits, fn, by = "species")

# save
# write_csv(traits, file = "data/scorp-traits.csv")


# Env. Variables ------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)

# read in community data
comm <- readxl::read_excel("data/scorp-comm.xlsx", sheet = 1, skip = 1)

# get and idea of the sampling extent
# plot(comm$lon, comm$lat)
# -60 -20 # lon
# -4 -12 # lat

# read in raster files
rr <- rast(list.files("../../../Documents/gis/paleoclim-current/2_5min/", pattern = "\\.tif$", full.names = T))

# soil raster (https://gaez.fao.org/pages/hwsd)
sr <- rast("/Users/stenio/Documents/gis/HWSD2_RASTER/HWSD2.bil") # home pc

# soil raster metadata (to recover the soil classes)
mm <- readxl::read_excel("/Users/stenio/Documents/gis/HWSD2_RASTER/HWSD2-SMU.xlsx", skip = 1)

# resample soil raster to bioclimatic rasters
sr <- resample(sr, rr$bio_1, method = "near")

# spatial points
pts <- st_as_sf(
  data.frame(
    lon = comm$lon,
    lat = comm$lat
  ),
  coords = c("lon", "lat"),
  crs = 4326
)

# must be TRUE
identical(crs(rr), crs(pts))

# extract bioclimatic variables
vals <- as_tibble(extract(rr, pts))

# extract soil variables
sv <- as_tibble(extract(sr, pts))

# merge soil and bioclimatic data
vals <- left_join(x = vals, y = sv, by = "ID")

# no missing values
vals[!complete.cases(vals), ]

# recover dominant soil class
mm |>
  filter(HWSD2_SMU_ID %in% vals$HWSD2) |>
  select(HWSD2_SMU_ID, WRB2) |>
  rename("HWSD2" = "HWSD2_SMU_ID") -> mm

# add to the data
left_join(x = vals, y = mm, by = "HWSD2") |>
  select(-c(ID, HWSD2)) |>
  rename("soil" = "WRB2") |>
  relocate("soil", .before = "bio_1") -> vals

# save
x <- cbind(comm[, c("site", "forest", "lon", "lat")], vals)
# write_csv(x, file = "data/site-envs.csv")


# Phylogeny -----------------------------------------------------------------------------------

library(tidyverse)
library(phytools)
library(ggtree)

# read in community matrix
comm <- readxl::read_excel("data/scorp-comm.xlsx", sheet = 1, skip = 1)

# read in tree
phy <- read.tree("data/scorp-complete-tree.tre")

# species not in the tree
setdiff(colnames(comm)[grepl("_", colnames(comm))], phy$tip.label)

# bind Bothriurus rabeca to the clade asper-rochai
nd <- getMRCA(phy = phy, tip = phy$tip.label[grepl("Bothriurus_asper|rochai", phy$tip.label)])
phy <- phytools::bind.tip(tree = phy, tip.label = "BOTH_Bothriurus_rabeca", where = nd)

# bind Tityus martinpaechi to clade serrulatus-kuryi-stigmurus
nd <- getMRCA(phy = phy, tip = phy$tip.label[grepl("serrulatus|kuryi|stigmurus", phy$tip.label)])
phy <- phytools::bind.tip(tree = phy, tip.label = "BUTH_Tityus_martinpaechi", where = nd)

# bind Ananteris otavianoi and Ananteris sp. caetes to the clade franckei-mauryi
nd <- getMRCA(phy = phy, tip = phy$tip.label[grepl("Ananteris_franckei|mauryi", phy$tip.label)])
phy <- phytools::bind.tip(tree = phy, tip.label = "BUTH_Ananteris_otavianoi", where = nd)
nd <- getMRCA(phy = phy, tip = phy$tip.label[grepl("Ananteris_franckei|mauryi", phy$tip.label)])
phy <- phytools::bind.tip(tree = phy, tip.label = "BUTH_Ananteris_sp_caetes", where = nd)

# adjust tree for potential rounding issues
phy <- phytools::force.ultrametric(phy, method = "extend")

# prune the tree
phy <- drop.tip(phy = phy, tip = setdiff(phy$tip.label, colnames(comm)[grepl("_", colnames(comm))]))

# check one more time (must be 0)
setdiff(colnames(comm)[grepl("_", colnames(comm))], phy$tip.label)
rm(nd)

# check the pruned tree
ggtree(phy) +
  geom_tiplab(size = 2.6) +
  xlim(0, 500)

# save
# write.tree(phy, file = "data/scorp-pruned-tree.tre")
