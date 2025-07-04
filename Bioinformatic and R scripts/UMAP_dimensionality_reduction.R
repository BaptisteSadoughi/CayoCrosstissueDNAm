# UMAP projection
# Useful tutorial
# UMAP https://cran.r-project.org/web/packages/umap/vignettes/umap.html
# Corrplot https://taiyun.github.io/corrplot/

rm(list = ls())

library_list <- c("umap","RColorBrewer","svglite","tidyverse","rlang", "corrplot")

lapply(library_list, require, character.only=TRUE)

# Percent methylation matrices of fully covered regions for each tissue
pmeth = readRDS("/scratch/sbaptis7/Cayotissue_CpG_coverage/combined_CpG_Regions/full_matrices/Regions_full_pmeth14.rds")

# Convert matrices to data frames for merging and create a new column with row names
list_of_dfs <- lapply(pmeth, function(x) {
  df <- as.data.frame(x)
  df$rownames <- row.names(df)
  return(df)
})

# Define function for joining two data frames
join_two_dfs <- function(df1, df2) {
  full_join(df1, df2, by = "rownames")
}

# Reduce list of data frames to a single data frame by joining
merged_df <- Reduce(join_two_dfs, list_of_dfs)

# Keep only rows with no missing values
pmeth <- merged_df[complete.cases(merged_df),]
rownames(pmeth) <- pmeth$rownames
pmeth <- pmeth %>% select(-rownames)
pmeth = as.matrix(pmeth)

# Load metadata
metadata_lid = read.table("/scratch/sbaptis7/Cayo_meth_metadata/metadata_final_lidpids_Nov24.txt", sep = "\t", header = TRUE) %>% filter(lid_pid != "LID_109490_PID_10416")

pmeth <- pmeth[,colnames(pmeth) %in% metadata_lid$lid_pid]

metadata_lid <- metadata_lid[match(colnames(pmeth),metadata_lid$lid_pid),]

identical(colnames(pmeth), metadata_lid$lid_pid)

#############################################
#     UMAP
#############################################

# Customized umap using Chiou et al. parameters https://www.nature.com/articles/s41593-022-01197-0
umap.custom <- umap.defaults
umap.custom$n_neighbors <- 50
umap.custom$min_dist <- 0.5

# Initialize list to store umap results and layout outputs
umap_results <- list()
umap_coords <- list()

set.seed(3500)
# Run UMAP with different configurations
for (config_name in c("defaults", "custom")) {
  umap_output <- umap(t(pmeth), config = get(paste0("umap.",config_name)))
  umap_results[[config_name]] <- umap_output  # Save umap results
  umap_coords[[config_name]] <- umap_output$layout  # Extract layout output
}

### PLOTTING

# Define the updated color order based on desired tissue arrangement
color_order <- c("adrenal", "heart", "kidney", "liver", "lung", "omental_at", 
                 "ovaries", "pituitary", "skeletal_muscle", "spleen", "testis", 
                 "thymus", "thyroid", "whole_blood")

# Sort tissue_plot and include the new levels
tissue_plot <- sort(c("whole_blood", "spleen", "omental_at", "heart", "kidney", 
                      "lung", "adrenal", "thymus", "thyroid", "pituitary", 
                      "liver", "skeletal_muscle"))

# Create the extended palette for the tissues
extended_palette <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(12), tissue_plot)

# Add the new levels and their respective colors
new_levels <- c("ovaries", "testis")
new_colors <- c("#008B8B", "#4682B4")
extended_palette <- c(extended_palette, setNames(new_colors, new_levels))

# Reorder the extended_palette to match the color_order
extended_palette <- extended_palette[color_order]

# Ensure metadata uses the updated levels
metadata_lid$grantparent_tissueType <- factor(metadata_lid$grantparent_tissueType, levels = color_order)

# Function to generate a discrete color palette
get_palette <- function(n) {
  # Check if n is 2, return custom palette
  if (n == 2) {
    return(c("grey","red"))
  }
  # For n>2 choose a suitable palette name and adjust the number of colors if necessary
  if (n <= 3) {
    return(brewer.pal(n, "Set1"))
  } else if (n <= 5) {
    return(brewer.pal(n, "Set2"))
  } else if (n<=12){
    return(brewer.pal(n, "Set3"))
  }
  else {
    return(viridis::viridis(n))
  }
}

### FUNCTION for customed UMAP
plot.umap = function(x, metadata, fill, shape = NULL,
                     mainTitle=TRUE,
                     main="A UMAP visualization of the dataset",
                     pad=0.1, cex.point=0.6, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=0.85, cex.axisTitle=1, cex.axisText=1,
                     legendfill.pos = "topleft", legendshape.pos="left", box = TRUE, color_order = color_order, custom_palette=NULL) {
  
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  
  # Adjust x and y limits with padding
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  
  if (!add) {
    par(mar=c(5, 5, 4, 2) + 0.1)  # Adjust the margins for the axes
    plot(xylim, xylim, type="n", frame=F, xlab="UMAP1", ylab="UMAP2", cex.lab=cex.axisTitle, cex.axis=cex.axisText) 
    if (mainTitle) {
      title(main="Your Main Title") #Add main title if main=TRUE
    }
    if (box) {
      box()  # Add a box around the plot
    }
  }
  
  # Determine whether the fill variable is discrete or continuous
  is_fill_discrete <- is.factor(metadata[[fill]]) || is.character(metadata[[fill]])
  
  # If the variable is continuous, we sort and make levels based on its values
  levels_fill <- if (is.numeric(metadata[[fill]])) {
    sort(unique(metadata[[fill]]))
  } else {
    unique(metadata[[fill]])
  }
  
  fill_factor <- factor(metadata[[fill]], levels = color_order)  # Reorder the factor according to color_order
  
  # If custom_palette is provided, we use it for both discrete and continuous variables
  if (!is.null(custom_palette)) {
    if (is_fill_discrete) {
      # For discrete data, we use the custom_palette (ensure it is long enough)
      if (length(custom_palette) < length(levels_fill)) {
        stop("Error: custom_palette does not have enough colors for all levels of the fill factor.")
      }
      colors <- custom_palette[1:length(levels_fill)]
    } else {
      # For continuous data, we create a gradient based on custom_palette
      colors <- colorRampPalette(custom_palette)(length(levels_fill))  # Color gradient for continuous values
    }
  } else {
    # If custom_palette is not provided, fall back to default behavior
    if (is_fill_discrete) {
      # Use a default palette for discrete variables
      palette <- get_palette(length(levels_fill))
      colors <- setNames(palette, levels(fill_factor))
      colors <- colors[color_order]
    } else {
      # Use a default continuous color palette for continuous variables
      colors <- colorRampPalette(brewer.pal(9, "Blues"))(length(levels_fill))  # Generate continuous colors
    }
  }
  
  # For continuous predictors, we map the actual data values to the color palette
  fill_as_int <- as.integer(fill_factor)  # This gives the integer level for each point
  
  # Handle continuous color mapping by scaling the actual data values
  if (!is_fill_discrete) {
    fill_values <- metadata[[fill]]
    color_scale <- scales::rescale(fill_values, to = c(1, length(colors)))
    point_fill_color <- colors[round(color_scale)]  # Map scaled values to color range
  } else {
    point_fill_color <- colors[fill_as_int]  # Use the discrete palette directly for discrete variables
  }
  
  # Handle shape if provided (max 5 levels)
  if (!is.null(shape)) {
    shape_levels <- unique(metadata[[shape]])
    if (length(shape_levels) > 5) {
      stop("Shape variable has more than 5 levels. Consider another approach.")
    }
    shape_pch <- 20 + (1:length(shape_levels))
    pch_vector <- sapply(metadata[[shape]], function(x) shape_pch[which(shape_levels == x)])
  } else {
    pch_vector <- rep(21, length(layout[, 1]))
  }
  
  # Plot the points on the UMAP
  points(layout[, 1], layout[, 2], col=NA, pch=pch_vector, bg=point_fill_color, cex=cex.point)
  
  # Add the legend for fill
  if (!add && !is.null(legendfill.pos)) {
    if (is_fill_discrete) {
      # Discrete legend
      labels.u <- levels(fill_factor)
      fills.u <- colors  # Colors for the discrete categories
    } else {
      # Continuous legend
      min_val <- min(metadata[[fill]])  # Minimum value of the continuous variable
      max_val <- max(metadata[[fill]])  # Maximum value of the continuous variable
      labels.u <- seq(min_val, max_val, length.out = 10)  # Create 5 labels
      labels.u <- round(labels.u)  # Round the labels to integers
      
      # Generate colors for the legend labels based on the continuous data
      rescaled_values <- scales::rescale(labels.u, to = c(1, length(colors)))
      fills.u <- colors[round(rescaled_values)]  # Map to colors and round to integer indices
    }
    
    legend.text <- as.character(labels.u)  # Convert labels to character for the legend
    legend(legendfill.pos, legend=legend.text, inset=0.03,
           col=NA, pt.cex = cex.point, pch=21, cex=cex.legend, pt.bg=fills.u, bty="n")
  }
  
  # Add the legend for shape
  if (!is.null(shape) && !is.null(legendshape.pos)) {
    shapes <- unique(metadata[[shape]])
    legend(legendshape.pos, legend = as.character(shapes), pch = shape_pch,
           bty = "n", pt.cex = cex.point, cex = cex.legend)
  }
  
  # Add main title
  if(mainTitle){
    mtext(side=3, main, cex=cex.main)
  }
}

## Generate plot
predictors <- c("grantparent_tissueType")

if(!dir.exists("/path/to/UMAP")) dir.create("/path/to/UMAP")

for(i in predictors){
  for (config_name in names(umap_coords)) {

    svglite::svglite(paste0("/path/to/UMAP/UMAP_",i,"_",config_name,".svg"))
    
    #if i is the tissues grab custom_palette otherwise set to NULL
    if(i == "grantparent_tissueType"){
      custom_palette=extended_palette
    } else {
        custom_palette = NULL
        }
    
    # Set number of plots per row and column
    par(mfrow = c(1,1))  # Replace 1 and 2 with desired number of rows and columns

    layout <- umap_coords[[config_name]]
    plot.umap(layout, metadata_lid, fill=i, shape=NULL, mainTitle = FALSE,
              cex.axisTitle = 1.7, cex.axisText = 1.7,cex.legend=0.01,legendfill.pos = NULL,
              main = NULL, box = FALSE, color_order = color_order, custom_palette = custom_palette)
    dev.off()
  }
}
