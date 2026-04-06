###########################################################################
# Population genetics reveals insights into Cryophytum} biogeography in South Africa
# van Steenderen, C.J.M., Sandenbergh, E., and Paterson, I.D. 2026
# Ecology and Evolution

# This script works through the climate analysis, focusing on WorldClim
# and soil variables, and linking that to the genetic SNP data
###########################################################################

source("setup.R")

#############################################################################
# CURRENT CLIMATE
#############################################################################

# Note: download these WorldClim and soil variables using geodata before proceeding
# due to size constraints, these files are not provided

pred_clim_current = terra::rast( list.files(
  here::here("climate_layers/current/wc2.1_2.5m/") ,
  full.names = TRUE,
  pattern = '.tif'
))  

# set the CRS (coordinate reference system) projection for the current climate layers
terra::crs(pred_clim_current) = "epsg:4326"
terra::crs(pred_clim_current, describe = T)

names(pred_clim_current) = paste0("bio_", sub(".*_", "", names(pred_clim_current)))

sample.info = readxl::read_excel("data/RADseq_sample_info.xlsx",
                                 sheet = 1) %>% janitor::clean_names()

str(sample.info)

sample.info$lat = as.numeric(sample.info$lat)
sample.info$long = as.numeric(sample.info$long)

#keep just SA

sample.info = sample.info %>%
  dplyr::filter(country == "South Africa")

# Convert to a spatial object for extraction
sample.pts = vect(sample.info, geom = c("long", "lat"), crs = "EPSG:4326")

# Extract climate values for each sample point
clim_values = terra::extract(pred_clim_current, sample.pts)

clim_values$site = sample.info$site
clim_values$sample_id = sample.info$sample_id
clim_values$country = sample.info$country
clim_values$splitstree = sample.info$splitstree

head(clim_values)

##########################
# Get soil data
#########################

###################################
# Download total nitrogen in soil 
###################################

# pred_soil = geodata::soil_world(
#   var = c("bdod", "cfvo", "clay", "nitrogen", "ocd", "phh2o", "sand", "silt", "soc"),
#   # Which variable do we want?
#   depth = 5,
#   # Soil depth (cm)
#   stat = "mean",
#   # Return mean values (could get variance, CI's, ect...)
#   path = here::here("soil_layers/")
# )

# load soil rasters
soil_rasters = terra::rast(list.files(
  here::here("soil_layers/soil_world/") ,
  full.names = TRUE,
  pattern = '.tif'
))

terra::crs(soil_rasters) = "epsg:4326"

# extract values for each site
soil_vals = terra::extract(soil_rasters, sample.pts)

environ.values = dplyr::bind_cols(clim_values, soil_vals)

environ.values = environ.values %>%
  dplyr::group_by(site) %>%
  dplyr::slice(1) %>%          # keep the first row per site
  dplyr::ungroup()

names(environ.values)

# PCA

# a-priori variable selection
selected_vars = c("bio_1", "bio_5", "bio_6", "bio_12",
                   "cfvo_0-5cm", "clay_0-5cm_mean", "nitrogen_0-5cm", "ocd_0-5cm",
                   "phh2o_0-5cm", "sand_0-5cm", "silt_0-5cm", "soc_0-5cm")

clim_pca_data = environ.values %>%
  dplyr::select(site, splitstree, dplyr::all_of(selected_vars)) %>% 
  na.omit() %>%
  filter(site != "MGSA")

head(clim_pca_data)

# check for collinearity

num_data = clim_pca_data %>%
  dplyr::select(where(is.numeric))

cor_mat = cor(num_data, use = "complete.obs", method = "pearson")

round(cor_mat, 2)

write.csv(cor_mat, "figures/corrmatrix.csv", row.names = F)

corplot = corrplot(cor_mat, method = "circle", type = "upper", tl.cex = 0.7, addCoef.col = "black")

png("figures/corrplot.png", width = 2000, height = 2000, res = 300)
corrplot(
  cor_mat,
  method = "circle",
  type = "upper",
  tl.cex = 0.7,
  addCoef.col = "black"
)
dev.off()

clim_pca_data = clim_pca_data %>%
  dplyr::select(-c(`ocd_0-5cm`, `silt_0-5cm`))

head(clim_pca_data)

# prepare the numeric matrix for PCA
clim_mat = clim_pca_data %>%
  dplyr::select(-c(site,splitstree)) %>%
  scale() 

# run PCA
pca_res = prcomp(clim_mat, center = TRUE, scale. = TRUE)
summary(pca_res)
pca_scores = pca_res$rotation

pca_scores

rownames(pca_scores) = c(
  "Mean Temp (Bio1)",
  "Max Temp (Bio5)",
  "Min Temp (Bio6)",
  "Annual Precip (Bio12)",
  "Coarse Fragments",
  "Clay %",
  "Soil N",
  #"Organic C",
  "Soil pH",
  "Sand %",
  # "Silt %",
  "Soil Organic C"
)

write.csv(pca_scores, "data/PCA.scores.csv")

# create a data frame with PCA scores and site info
pca_df = as.data.frame(pca_res$x) %>%
  mutate(site = clim_pca_data$site,
         splitstree = clim_pca_data$splitstree) %>%
  dplyr::select(PC1, PC2, PC3, site, splitstree)

pca_df

# plot PCA coloured by site
PCA.plot = ggplot(pca_df, aes(x = PC1, y = PC2, fill = splitstree)) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  labs(#title = "c) PCA",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"),
    fill = "SplitsTree clade") +
  theme(legend.position = "right") +
  theme_classic() +
  ggrepel::geom_text_repel(aes(label = site),
                           size = 3,         # text size
                           max.overlaps = 20, # adjust as needed
                           box.padding = 0.3, # distance from points
                           segment.color = 'grey50') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +  # horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")    # vertical line

PCA.plot

ggsave("figures/PCA.environs.png", PCA.plot, width = 7, height = 6, dpi = 350)

custom_var_labels = c(
  "Mean Temp (Bio1)",
  "Max Temp (Bio5)",
  "Min Temp (Bio6)",
  "Annual Precip (Bio12)",
  "Coarse Fragments",
  "Clay %",
  "Soil N",
  # "Organic C",
  "Soil pH",
  "Sand %",
  #  "Silt %",
  "Soil Organic C"
)

pca_res$rotation

rownames(pca_res$rotation) = custom_var_labels

biplot.pca = fviz_pca_biplot(
  pca_res,
  geom.ind = "point",
  palette = c("green", "pink", "cyan", "darkblue", "red", "gold"),
  fill.ind = pca_df$splitstree,
  pointshape = 21, 
  pointsize = 4,
  mean.point = FALSE,
  repel = TRUE,
  label = "var",  
  addEllipses = TRUE, 
  ellipse.type = "convex",
  label.var = list(labels = custom_var_labels),
  col.var = "grey50",
  legend.title = list(fill = "SplitsTree clade")
) +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  scale_color_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  theme_classic()

biplot.pca

# ggsave("figures/PCA.biplot.png", biplot.pca, 
#        width = 8, height = 6, dpi = 350)
# ggsave("figures/PCA.biplot.svg", biplot.pca, 
#        width = 8, height = 6, dpi = 350)

PCA.plot.2 = ggplot(pca_df, aes(x = PC1, y = PC3, fill = splitstree)) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  labs(
    #title = "c) PCA",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC3 (", round(summary(pca_res)$importance[2, 3] * 100, 1), "%)"),
    fill = "SplitsTree clade"
  ) +
  theme(legend.position = "right") +
  theme_classic() +
  ggrepel::geom_text_repel(
    aes(label = site),
    size = 3,           # text size
    max.overlaps = 20,  # adjust as needed
    box.padding = 0.3,  # distance from points
    segment.color = 'grey50'
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +  # horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")    # vertical line

PCA.plot.2

PCA.plot.3 = ggplot(pca_df, aes(x = PC2, y = PC3, fill = splitstree)) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  labs(
    #title = "c) PCA",
    x = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
    y = paste0("PC3 (", round(summary(pca_res)$importance[2, 3] * 100, 1), "%)"),
    fill = "SplitsTree clade"
  ) +
  theme(legend.position = "right") +
  theme_classic() +
  ggrepel::geom_text_repel(
    aes(label = site),
    size = 3,          
    max.overlaps = 20, 
    box.padding = 0.3, 
    segment.color = 'grey50'
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")    

PCA.plot.3

loadings_pc1 = pca_res$rotation[,1]  
loadings_pc1 = sort(abs(loadings_pc1), decreasing = TRUE) %>%  # largest contributors first
  as.data.frame()
colnames(loadings_pc1) = c("loading")
loadings_pc1$variable = rownames(loadings_pc1)
loadings_pc1$variable = factor(loadings_pc1$variable, levels = loadings_pc1$variable)
loadings_pc1$pc = "PC1"

loadings_pc1 = loadings_pc1 %>%
  dplyr::mutate(variable_name = recode(variable,
                                       "ocd_0-5cm"        = "Organic C",
                                       "soc_0-5cm"        = "Soil Organic C",
                                       "bio_5"            = "Max Temp (Bio5)",
                                       "phh2o_0-5cm"      = "Soil pH",
                                       "bio_12"           = "Annual Precip (Bio12)",
                                       "bio_6"            = "Min Temp (Bio6)",
                                       "nitrogen_0-5cm"   = "Soil N",
                                       "clay_0-5cm_mean"  = "Clay %",
                                       "bio_1"            = "Mean Temp (Bio1)",
                                       "sand_0-5cm"       = "Sand %",
                                       "silt_0-5cm"       = "Silt %",
                                       "cfvo_0-5cm"       = "Coarse Fragments"
  ))

loadings_pc2 = pca_res$rotation[,2]  
loadings_pc2 = sort(abs(loadings_pc2), decreasing = TRUE) %>%  # largest contributors first
  as.data.frame()
colnames(loadings_pc2) = c("loading")
loadings_pc2$variable = rownames(loadings_pc2)
loadings_pc2$variable = factor(loadings_pc2$variable, levels = loadings_pc2$variable)
loadings_pc2$pc = "PC2"

loadings_pc2 = loadings_pc2 %>%
  dplyr::mutate(variable_name = recode(variable,
                                       "sand_0-5cm"       = "Sand %",
                                       "silt_0-5cm"       = "Silt %",
                                       "nitrogen_0-5cm"   = "Soil N",
                                       "clay_0-5cm_mean"  = "Clay %",
                                       "cfvo_0-5cm"       = "Coarse Fragments",
                                       "bio_12"           = "Annual Precip (Bio12)",
                                       "phh2o_0-5cm"      = "Soil pH",
                                       "bio_1"            = "Mean Temp (Bio1)",
                                       "bio_6"            = "Min Temp (Bio6)",
                                       "soc_0-5cm"        = "Soil Organic C",
                                       "ocd_0-5cm"        = "Organic C",
                                       "bio_5"            = "Max Temp (Bio5)"
  ))

head(loadings_pc1)
head(loadings_pc2)

loadings.df = dplyr::bind_rows(loadings_pc1, loadings_pc2)

# Plot
PC1.plot = ggplot(loadings_pc1, aes(x = variable_name, y = loading, fill = loading)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "a) PC1",
       x = "Abiotic Variable",
       y = "Absolute PCA Loading") +
  theme_classic() +
  theme(legend.position = "none") 

PC2.plot = ggplot(loadings_pc2, aes(x = variable_name, y = loading, fill = loading)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "b) PC2",
       x = "Abiotic Variable",
       y = "Absolute PCA Loading") +
  theme_classic() +
  theme(legend.position = "none")

loading.plots = gridExtra::grid.arrange(PC1.plot, PC2.plot, ncol = 2)

ggsave("figures/PCA.loadings.png", loading.plots, width = 9, height = 3, dpi = 350)


##################################
# SNP pca
##################################

all.snps = vcfR::read.vcfR("data/populations.snps.filtered_mac3.vcf")
all.snps_genlight = vcfR::vcfR2genlight(all.snps)
genlight_sample_names = as.data.frame(all.snps_genlight@ind.names)
colnames(genlight_sample_names) = "sample_id"

samples_to_remove = indNames(all.snps_genlight)[
  grepl("^(USA|MX|EG|CYP|I|MO|CI)", indNames(all.snps_genlight))
]

# Subset the genlight object to exclude them
all.snps_genlight.filtered = all.snps_genlight[!indNames(all.snps_genlight) %in% samples_to_remove, ]

# set populations
all.snps_genlight.filtered@pop = as.factor(all.snps_genlight.filtered@ind.names)

sample.info.filtered = sample.info %>%
  dplyr::filter(sample_id %in% all.snps_genlight.filtered@ind.names)

sample.info.filtered = sample.info.filtered[match(all.snps_genlight.filtered@ind.names, 
                                                  sample.info.filtered$sample_id), ]

genlight_sample_names$sample_id = trimws(genlight_sample_names$sample_id)

filtered_info = sample.info %>%
  dplyr::filter(sample_id %in% genlight_sample_names$sample_id)

# Reorder to match genlight_sample_names
filtered_info = filtered_info %>%
  dplyr::slice(match(genlight_sample_names$sample_id, sample_id))

str(filtered_info)

PCA = adegenet::glPca(all.snps_genlight, nf = 3)

scores = as.data.frame(PCA$scores)  # PCA coordinates (individuals × axes)

head(scores)

# Remove rows where the rowname contains "USA" or "MX"
scores = scores[
  (grepl("SA", rownames(scores)) | grepl("OUSA", rownames(scores))) &
    !grepl("^USA", rownames(scores)) &  # Exclude names starting with USA
    !grepl("MX", rownames(scores)),     # Exclude anything with MX
]

scores$sample_id = rownames(scores)
head(scores)

scores = scores %>%
  dplyr::left_join(sample.info %>% 
                     dplyr::select(sample_id, site, splitstree),
                   by = "sample_id")

head(scores)
scores

scores = scores %>%
  dplyr::distinct(site, .keep_all = TRUE)


# Environmental PCA scores
pca_df
pca_df$type = "envs"
head(pca_df)

# SNP PCA scores
scores
scores$type = "snps"
head(scores)

scores = scores %>%
  dplyr::select(-c(sample_id, PC3))

# Select relevant columns
env_site = pca_df %>%
  dplyr::filter(type == "envs") %>%
  dplyr::select(site, splitstree, PC1_env = PC1, PC2_env = PC2)

snp_site = scores %>%
  dplyr::filter(type == "snps") %>%
  dplyr::select(site, splitstree, PC1_snp = PC1, PC2_snp = PC2)

# Merge by site
combined = dplyr::inner_join(env_site, snp_site, by = "site")

combined

# guerichianum only
combined_guerich = combined %>% dplyr::filter(splitstree.x != "crystallinum", splitstree.x != "hybrid")  %>%
  dplyr::left_join(
    sample.info %>%
      #dplyr::mutate(site = substr(sample_id, 1, 4)) %>%  # extract site code (e.g., DWSA1 → DWSA)
      dplyr::distinct(site, .keep_all = TRUE) %>%        # keep only the first occurrence of each site
      dplyr::select(site, lat, long),                    # keep only relevant columns
    by = "site"
  )

# crystallinum only
combined_cryst = combined %>% dplyr::filter(splitstree.x == "crystallinum")  %>%
  dplyr::left_join(
    sample.info %>%
      #dplyr::mutate(site = substr(sample_id, 1, 4)) %>%  # extract site code (e.g., DWSA1 → DWSA)
      dplyr::distinct(site, .keep_all = TRUE) %>%        # keep only the first occurrence of each site
      dplyr::select(site, lat, long),                    # keep only relevant columns
    by = "site"
  )

combined = combined %>%
  dplyr::left_join(
    sample.info %>%
      #dplyr::mutate(site = substr(sample_id, 1, 4)) %>%  # extract site code (e.g., DWSA1 → DWSA)
      dplyr::distinct(site, .keep_all = TRUE) %>%        # keep only the first occurrence of each site
      dplyr::select(site, lat, long),                    # keep only relevant columns
    by = "site"
  )

# Check
head(combined)

cor.test(combined$PC1_env, combined$PC1_snp)
cor.test(combined$PC2_env, combined$PC2_snp)

lm_PC1 = lm(PC1_snp ~ PC1_env, data = combined)
lm_PC2 = lm(PC2_snp ~ PC2_env, data = combined)

# Extract coefficients and R-sq
coef_PC1 = coef(lm_PC1)
r2_PC1 = summary(lm_PC1)$r.squared

coef_PC2 = coef(lm_PC2)
r2_PC2 = summary(lm_PC2)$r.squared

eq_PC1 = paste0("y = ", round(coef_PC1[2], 2), "x + ", round(coef_PC1[1], 2),
                 "\nR² = ", round(r2_PC1, 2))
eq_PC2 = paste0("y = ", round(coef_PC2[2], 2), "x + ", round(coef_PC2[1], 2),
                 "\nR² = ", round(r2_PC2, 2))


PC1.corplot = ggplot(combined, aes(x = PC1_env, y = PC1_snp, fill = splitstree.x)) +
  geom_smooth(method = "lm", se = TRUE, colour = "royalblue", fill = "lightblue") +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  ggrepel::geom_text_repel(aes(label = site), size = 2) +
  theme_classic() +
  labs(
    title = "a)",
    subtitle = eq_PC1,
    x = "WorldClim + Soil PC1",
    y = "SNPs PC1",
    fill = "SplitsTree clade"
  ) +
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +  # horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") ; PC1.corplot

PC2.corplot = ggplot(combined, aes(x = PC2_env, y = PC2_snp, fill = splitstree.x)) +
  geom_smooth(method = "lm", se = TRUE, colour = "royalblue", fill = "lightblue") +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("green", "pink", "cyan", "darkblue", "red", "gold")) +
  ggrepel::geom_text_repel(aes(label = site), size = 2) +
  theme_classic() +
  labs(
    title = "b)",
    subtitle = eq_PC2,
    x = "WorldClim + Soil PC2",
    y = "SNPs PC2",
    fill = "SplitsTree clade"
  ) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +  # horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") ; PC2.corplot   # vertical line

PC.corplots = gridExtra::grid.arrange(PC1.corplot, PC2.corplot, ncol = 2)
ggsave("figures/PCA.correlations.guerich.png", PC.corplots, 
       width = 10, height = 4, dpi = 350)

#################################################################
# Mantel tests
#################################################################

# get GPS points

head(combined_cryst)

# remove YZSA to see - the points at the top of the plot that look different are all
# comparisons to that site

# NB: change the input here (combined, combined_guerich, or combined_cryst)

combined_guerich = combined_guerich %>%
  dplyr::filter(site != "YZSA")

coords = combined_guerich[, c("long", "lat")]
geo_dist = geosphere::distm(coords, fun = geosphere::distHaversine) / 1000  # convert meters to km
rownames(geo_dist) = combined_guerich$site
colnames(geo_dist) = combined_guerich$site

mantel(
  dist(combined_guerich[, c("PC1_env", "PC2_env")]),
  dist(combined_guerich[, c("PC1_snp", "PC2_snp")])
)

# Compute pairwise distances
env_dist = as.matrix(dist(combined_guerich[, c("PC1_env", "PC2_env")]))
snp_dist = as.matrix(dist(combined_guerich[, c("PC1_snp", "PC2_snp")]))

# partial mantel, accounting for distance between sites
mantel.partial(
  as.dist(snp_dist),
  as.dist(env_dist),
  as.dist(geo_dist),
  method = "pearson",
  permutations = 999
)

# Convert matrices to long format
env_long = as.data.frame(as.table(env_dist))
snp_long = as.data.frame(as.table(snp_dist))

# Merge the two
dist_df = env_long %>%
  rename(site1 = Var1, site2 = Var2, env_dist = Freq) %>%
  left_join(snp_long %>% rename(site1 = Var1, site2 = Var2, snp_dist = Freq),
            by = c("site1", "site2")) %>%
  # Remove self-comparisons (distance = 0)
  filter(site1 != site2)

summary(lm(data = dist_df, snp_dist ~ env_dist))

mantel.plot = ggplot(dist_df, aes(x = env_dist, y = snp_dist)) +
  geom_smooth(method = "lm", color = "grey60", fill = "lightgrey", se = TRUE, lwd = 0.6) +
  geom_point(alpha = 0.6, color = "darkblue") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +
  theme_classic() +
  labs(
    x = "Environmental distance (PC1+PC2)",
    y = "Genetic distance (PC1+PC2)"
  ) ;mantel.plot

#ggsave("figures/mantel.guerich.png", mantel.plot, width = 6, height = 4, dpi = 350)

#####################################################################################
# END
#####################################################################################