################################################################
# PLOT STRUCTURE/fastSTRUCTURE output
################################################################

# inputs:
# output file from fastSTRUCTURE (.meanQ)
# sample info data frame, with two columns: id and pop

#################################################################

structure.plot = function(strucfile, sampleinfo, kval="N"){
  
  f <- readLines(strucfile, warn = FALSE)
  
  qmat <- read.delim(
    strucfile, 
    header = FALSE,
    sep = "",
    strip.white = TRUE
  )
  
  names(qmat) <- c(paste0("pop_",1:(ncol(qmat))))
  
  qmat$name = sampleinfo$id
  qmat$population = sampleinfo$pop
  
  qmat$population <- as.factor(qmat$population)
  qmat <- qmat %>% 
    pivot_longer(c(-name, -population), names_to = "group", values_to = "probability")
  
  my_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  my_pal = c("#66A61E", "yellow", "lightblue", "#D95F02", "darkblue", "#666666", "#7570B3", "#A6761D",
             "#E7298A", "red" )
  

  struc.plot = qmat %>% 
    ggplot(aes(x = name, y = probability, fill = group)) +
    geom_bar(stat = "identity", width = 1.0) +
    theme_bw() +
    labs(title = paste0("K = ", kval)) +
    ylab("Probability") +
    xlab("") +
    #guides(fill=guide_legend(title="Membership")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off") +
    scale_fill_manual(values=c(my_pal)) +
    facet_grid(~population, scales = "free_x", space = "free" )
  
  return(struc.plot)
  
}#structure.plot



hapmap.global = function(struc.file, sample.info, gps.data, piesize = 0.85){
  
  f <- readLines(struc.file, warn = FALSE)
  
  qmat <- read.delim(
    struc.file, 
    header = FALSE,
    sep = "",
    strip.white = TRUE
  )
  
  names(qmat) <- c(paste0("pop_",1:(ncol(qmat))))
  
  qmat$name = sample.info$id
  qmat$population = sample.info$pop
  
  qmat$population <- as.factor(qmat$population)
  qmat <- qmat %>% 
    pivot_longer(c(-name, -population), names_to = "group", values_to = "probability")
  
  #my_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  my_pal = c("#66A61E", "yellow", "lightblue", "#D95F02", "darkblue", "#666666", "#7570B3", "#A6761D",
             "#E7298A", "red" )
  
  qmat_merged <- qmat %>%
    group_by(name) %>%
    mutate(lat = gps.data$lat[match(name, gps.data$name)],
           lon = gps.data$lon[match(name, gps.data$name)]) %>%
    ungroup()
  
  qmat_merged$lat = as.numeric(qmat_merged$lat)
  qmat_merged$lon = as.numeric(qmat_merged$lon)
  
  # this is a summary for each individual sample
  qmat_wide <- qmat_merged %>%
    tidyr::pivot_wider(names_from = group, values_from = probability)
  
  # this averages all the data per unique group, across all samples for that group
  qmat_averaged <- qmat_wide %>%
    group_by(population) %>%
    summarize(
      lat = first(lat),      # Keep the first value of lat
      lon = first(lon),      # Keep the first value of lon
      across(starts_with("pop"), mean, na.rm = TRUE)  # Average pop_* columns
    ) %>%
    dplyr::rename(name = population) 
  
  qmat_averaged <- qmat_averaged %>%
    mutate(
      lat = ifelse(is.na(lat), 35.0, lat),  # Replace NA lat with Cyprus lat
      lon = ifelse(is.na(lon), 33.0, lon)   # Replace NA lon with Cyprus lon
    )
  
  world_map <- rnaturalearth::ne_countries(
    scale = "medium", 
    returnclass = "sf"
  ) 
  
  pies = ggplot() +
    
    geom_sf(data = world_map, alpha = 0.2, fill = "grey80") +
    
    scatterpie::geom_scatterpie(
      #data = na.omit(qmat_wide), # plot pies for each individual sample
      data = na.omit(qmat_averaged), # plot pies representing the average per group
      aes(x = lon, y = lat, group = name),
      #aes(x = lon, y = lat, group = name),
      cols = pop_cols <- grep("^pop_", colnames(qmat_wide), value = TRUE),
      #alpha = 0.3,
      #color = NA, # Remove borders for pie charts
      color = "black", # Remove borders for pie charts
      linewidth = 0.2,
      pie_scale = piesize # Adjust pie size if necessary
    ) +
    
    # uncomment if you want labels on the pies
    # geom_text(
    #   data = na.omit(qmat_averaged),
    #   aes(x = lon, y = lat, label = name), # Assuming 'name' is the group label
    #   color = "black",
    #   fontface = "bold",
    #   nudge_y = 0.1 # Adjust this value if needed
    # ) +
    
    scale_fill_manual(values = my_pal) +
    
    coord_sf(
      ylim = c(95,-60),
      crs = 4326,
      expand = FALSE
    ) + 
    
    xlab("Longitude") + 
    ylab("Latitude") +
    theme_classic() +
    theme(legend.position = "none")
  
  return(pies)
} #hapmap.global
