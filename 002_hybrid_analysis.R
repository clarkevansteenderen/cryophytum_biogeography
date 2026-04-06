###########################################################################
# Population genetics reveals insights into Cryophytum} biogeography in South Africa
# van Steenderen, C.J.M., Sandenbergh, E., and Paterson, I.D. 2026
# Ecology and Evolution

# This script works through futher downstream data processing
###########################################################################

source("setup.R")
BiocManager::install("SNPRelate")

all.snps = vcfR::read.vcfR("data/populations.snps.filtered_mac3.vcf")
all.snps_genlight = vcfR::vcfR2genlight(all.snps)
genlight_sample_names = as.data.frame(all.snps_genlight@ind.names)
colnames(genlight_sample_names) = "sample_id"

sample.info = readxl::read_excel("data/RADseq_sample_info.xlsx",
                                 sheet = 1) %>% janitor::clean_names()

samples_to_remove = adegenet::indNames(all.snps_genlight)[
  grepl("^(USA|MX|EG|CYP|I|MO|CI)", indNames(all.snps_genlight))
]

# Subset the genlight object to exclude them
all.snps_genlight.filtered = all.snps_genlight[!indNames(all.snps_genlight) %in% samples_to_remove, ]

# set populations
all.snps_genlight.filtered@pop = as.factor(all.snps_genlight.filtered@ind.names)

# Verify
table(pop(all.snps_genlight.filtered))

# save a version for SplitsTree
sample.div = StAMPP::stamppNeisD(all.snps_genlight.filtered, pop = FALSE)

#export for splitstree
StAMPP::stamppPhylip(distance.mat=sample.div, 
                     file="data/splitstree.txt")

# save for faststructure
dartR::gl2faststructure(x = all.snps_genlight.filtered, outpath=getwd(),
                        outfile = "data/faststructure_input.str")

all.snps_genlight.filtered@ind.names
head(sample.info)

# subset the sample info so that it has the same samples as the filtered genlight

sample.info.filtered = sample.info %>%
  dplyr::filter(sample_id %in% all.snps_genlight.filtered@ind.names)

sample.info.filtered = sample.info.filtered[match(all.snps_genlight.filtered@ind.names, 
                                            sample.info.filtered$sample_id), ]


#write.csv(sample.info.filtered, "data/sample_info.csv", row.names = T)

# assign the two different pop structures: by clade (as seen in splitstree) and by site

all.snps.clade = all.snps_genlight.filtered
all.snps.clade@pop = as.factor(sample.info.filtered$splitstree)
all.snps.clade = dartR::gl2gi(all.snps.clade) 
  
all.snps.site = all.snps_genlight.filtered
all.snps.site@pop = as.factor(sample.info.filtered$site)
all.snps.site = dartR::gl2gi(all.snps.site) 

# FST values -> only works on genind objects
clade.FST = hierfstat::genet.dist(all.snps.clade, method = "WC84")
clade.FST = as.data.frame(as.matrix(clade.FST)) %>%
  round(., digits = 2)

clade.FST.boots = hierfstat::boot.ppfst(all.snps.clade, nboot = 10000)
clade.FST.boots$ll = round(clade.FST.boots$ll, 2)
clade.FST.boots$ul = round(clade.FST.boots$ul, 2)

# Combine them into a single character matrix (e.g. "0.06/0.09")
ci_combined_clades = matrix(NA, nrow = nrow(clade.FST.boots$ll), ncol = ncol(clade.FST.boots$ll),
                      dimnames = dimnames(clade.FST.boots$ll))

for (i in 1:nrow(clade.FST.boots$ll)) {
  for (j in 1:ncol(clade.FST.boots$ll)) {
    if (!is.na(clade.FST.boots$ll[i, j]) && !is.na(clade.FST.boots$ul[i, j])) {
      ci_combined_clades[i, j] = paste0(clade.FST.boots$ll[i, j], "/", clade.FST.boots$ul[i, j])
    }
  }
}

ci_combined_clades

# Convert clade.FST to character so we can combine numeric and string values
clade.FST.labeled = as.matrix(clade.FST)
clade.FST.labeled = format(round(clade.FST.labeled, 2), nsmall = 2)  # keep formatting

# Replace the upper triangle with CI values from ci_combined_clades
upper_idx = upper.tri(clade.FST.labeled)

clade.FST.labeled[upper_idx] = ci_combined_clades[upper_idx]

clade.FST.labeled
clade.FST.labeled = as.data.frame(clade.FST.labeled)

#write.csv(clade.FST.labeled, "publication_hybrid/fst_clades.csv", row.names = T)

################################################
# Pop stats
################################################

###################################
# STANDARD PCA
###################################

sample.info = readxl::read_excel("data/RADseq_sample_info.xlsx",
                                 sheet = 1) %>% janitor::clean_names()

head(sample.info)
names(sample.info)
str(sample.info)

sample.info$sample_id = trimws(sample.info$sample_id)

genlight_sample_names$sample_id = trimws(genlight_sample_names$sample_id)

filtered_info = sample.info %>%
  dplyr::filter(sample_id %in% genlight_sample_names$sample_id)

# Reorder to match genlight_sample_names
filtered_info = filtered_info %>%
  dplyr::slice(match(genlight_sample_names$sample_id, sample_id))

str(filtered_info)

PCA = adegenet::glPca(all.snps_genlight, nf = 3)

# Extract PCA scores
scores = as.data.frame(PCA$scores)  # PCA coordinates (individuals × axes)
scores$group = filtered_info$biogeography

head(scores)

# Remove rows where the rowname contains "USA" or "MX"
scores = scores[
  (grepl("SA", rownames(scores)) | grepl("OUSA", rownames(scores))) &
    !grepl("^USA", rownames(scores)) &  # Exclude names starting with USA
    !grepl("MX", rownames(scores)),     # Exclude anything with MX
]

####################################
# PCA
####################################

as.factor(scores$group)

pca.ggplot = ggplot(scores, aes(x = PC1, y = PC2, fill = group)) +
  
  ggrepel::geom_text_repel(aes(label = rownames(scores)), 
                           vjust = -1.5, size = 2.5, colour = "black",
                           segment.color = "grey80",
                           max.overlaps = 50) +
  geom_point(alpha = 1, size = 3, shape = 21) +
  labs(
    #title = "PCA of RADseq data - ref genome assembly",
    x = paste0("PC1 (", round(PCA$eig[1] / sum(PCA$eig) * 100, 1), "%)"),
    y = paste0("PC2 (", round(PCA$eig[2] / sum(PCA$eig) * 100, 1), "%)")
  ) +
  scale_fill_manual(values = c(
    "blue",
    "lightblue",
    "red",
    "gold",
    "forestgreen",
    "green"
  )) +
  theme_classic() +
  theme(legend.position = "right") +
  #ggtitle("d)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

pca.ggplot

# ggsave("figures/pca.png", pca.ggplot, width = 12, dpi=400,
#        height = 7, units = "in")
# ggsave("figures/pca.svg", pca.ggplot, width = 11, dpi=400,
#        height = 7, units = "in")

#################################################################
# MAP OF SA
#################################################################

gps.points = readxl::read_excel("data/RADseq_sample_info.xlsx",
                                sheet = 1) %>%
  janitor::clean_names() %>%
  dplyr::select(sample_id, lat, long, biogeography, splitstree) 

gps.points$lat = as.numeric(gps.points$lat)
gps.points$long = as.numeric(gps.points$long)

southafrica_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name %in% c("South Africa", "Lesotho", "eSwatini"))

provincial_borders = rnaturalearth::ne_states(country = "South Africa", returnclass = "sf")

# Plot GPS points on world map to check our locality data is correct 
distr_map = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = southafrica_ext, alpha = 0.5, fill = "grey80") +
  # Add GPS points 
  
  # ggrepel::geom_text_repel(
  #   data = gps.points,
  #   aes(x = long, y = lat, label = sample_id),  # replace `id` with your sample ID column
  #   size = 2,
  #   color = "grey15",
  #   max.overlaps = 100,   # adjust as needed
  #   box.padding = 0.25,
  #   point.padding = 0.1,
  #   segment.color = "grey85"
  # ) +
  geom_sf(data = provincial_borders, color = "black") +
  geom_point(
    data = gps.points, 
    aes(
      x = long, 
      y = lat, 
      fill = I(ifelse(biogeography == "crystallinum", "black", "white"))
    ), 
    size = 2,
    shape = 21,
    color = "black",
    alpha = 0.75
  ) +
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    xlim = c(16, 28),      
    ylim = c(-35, -28),
    expand = FALSE
  ) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  #ggtitle(expression(italic("Orachrysops niobe") * ", South Africa")) +
  ggspatial::annotation_scale(
    location = "br",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  theme_classic() +
  theme(legend.position = "top")

distr_map

# ggsave("figures/distr_map.png", distr_map, width = 10, dpi=400,
#        height = 5, units = "in")
# ggsave("figures/distr_map_materials_provs.svg", distr_map, width = 7, dpi=400,
#        height = 7, units = "in")

####################################################################################
# ELEVATION MAP
####################################################################################

elevation = geodata::elevation_global(res = 2.5, path = tempdir())

SA_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name %in% c("South Africa", "Lesotho", "eSwatini"))

elevation_SA = raster::mask(elevation, SA_ext)
elevation_SA_df = as.data.frame(elevation_SA, xy = TRUE)
colnames(elevation_SA_df) = c("lon", "lat", "altitude")
head(elevation_SA_df)

SA_elevation_map = ggplot() +
  geom_raster(data = elevation_SA_df, aes(x = lon, y = lat, fill = altitude)) +
  theme_classic() +
  geom_point(
    data = gps.points, 
    aes(x = long, y = lat), 
    size = 2,
    shape = 21,
    color = "black",
    fill = "white",
    alpha = 0.75
  ) +
  #scale_fill_viridis_c(name = "Altitude (m)") +
  scale_fill_gradientn(
    colours = MetBrewer::met.brewer("Demuth", type = "continuous"),
    name = "Altitude (m)"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    xlim = c(16, 28),      
    ylim = c(-35, -28),
    expand = FALSE
  ) +
  ggspatial::annotation_scale(
    location = "bl",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  )

SA_elevation_map

# ggsave("figures/elev_map.png", SA_elevation_map, width = 10, dpi=400,
#        height = 5, units = "in")
# ggsave("figures/elev_map.svg", SA_elevation_map, width = 7, dpi=400,
#        height = 7, units = "in")

#####################################################################################
# create STRUCTURE files

dartR::gl2faststructure(x = all.snps_genlight, outpath=getwd(),
                        outfile = "data/faststructure_input.str")

################################################################################

# faststructure results

source("structure.plot.R")

K2struc.file = "faststructure/faststructure.2.meanQ"
K3struc.file = "faststructure/faststructure.3.meanQ"
K4struc.file = "faststructure/faststructure.4.meanQ"
K5struc.file = "faststructure/faststructure.5.meanQ"
K6struc.file = "faststructure/faststructure.6.meanQ"

faststruc.popinfo = sample.info.filtered %>%
  dplyr::select(sample_id, splitstree) %>%
  dplyr::rename(id = sample_id,
                pop = splitstree)

faststruc.popinfo

# faststruc.popinfo = faststruc.popinfo %>%
#   dplyr::filter(pop != "crystallinum")

head(faststruc.popinfo)

K2.nocryst = structure.plot(K2struc.file, faststruc.popinfo, kval = "2"); K2.nocryst

K3.nocryst = structure.plot(K3struc.file, faststruc.popinfo, kval = "3"); K3.nocryst

K4.nocryst = structure.plot(K4struc.file, faststruc.popinfo, kval = "4"); K4.nocryst

K5.nocryst = structure.plot(K5struc.file, faststruc.popinfo, kval = "5"); K5.nocryst
K6 = structure.plot(K6struc.file, faststruc.popinfo, kval = "6"); K6

strucplots = gridExtra::grid.arrange(K2.nocryst, K3.nocryst, K4.nocryst, K5.nocryst,
                                     ncol = 1)

# ggsave("figures/structureplot.png", 
#        strucplots, width = 11, height = 8, units = "in")
# ggsave("figures/structureplot.svg", 
#        strucplots, width = 11, height = 8, units = "in")


# PIE CHARTS

faststruc.popinfo
f = readLines(K5struc.file, warn = FALSE)

qmat = read.delim(
  K5struc.file, 
  header = FALSE,
  sep = "",
  strip.white = TRUE
)

names(qmat) = c(paste0("pop_",1:(ncol(qmat))))

qmat$name = faststruc.popinfo$id
qmat$population = faststruc.popinfo$pop

qmat$population = as.factor(qmat$population)
qmat = qmat %>% 
  pivot_longer(c(-name, -population), names_to = "group", values_to = "probability")

#my_pal = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
my_pal = c("#66A61E", "yellow", "lightblue", "#D95F02", "darkblue", "#666666", "#7570B3", "#A6761D",
           "#E7298A", "red" )

qmat

qmat_wide = qmat %>%
  tidyr::pivot_wider(names_from = group, values_from = probability)

# this averages all the data per unique group, across all samples for that group
qmat_averaged = qmat_wide %>%
  group_by(population) %>%
  summarize(
    across(starts_with("pop"), mean, na.rm = TRUE)  # Average pop_* columns
  ) %>%
  dplyr::rename(name = population) 


# reshape data to long format
qmat_long = qmat_averaged %>%
  pivot_longer(cols = starts_with("pop"), names_to = "population", values_to = "proportion")

# plot pie charts, one per clade
K5.pies = ggplot(qmat_long, aes(x = "", y = proportion, fill = population)) +
  geom_col(width = 1, color = "white") +
  scale_fill_manual(values = my_pal) +
  coord_polar(theta = "y") +
  facet_wrap(~name) +
  theme_void() + 
  theme(legend.position = "bottom") +
  labs(fill = "Population")

  
ggsave("figures/struc.pies.svg", K5.pies, width = 8, height = 6, units = "in")

#####################################################################################
# END
#####################################################################################