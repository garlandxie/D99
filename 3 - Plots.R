# RLQ Analysis - Xie et al. 2018 
# Script is used to construct plots for RLQ analysis
# Code developed by Garland Xie

#### Libraries ####

# Load up the easypackages 
# checks to see if you have these packages. If not, install them
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# Load packages
packages("ggplot2", "ggrepel", "here", "stringr", "ade4")

# Load libraries
# ggplot2: data visualization
libraries("ggplot2", "ggrepel", "here", "stringr", "ade4")

#### Import #####
rlq_250 <- readRDS(here("Objects", "D99_rlq_250.rds"))
rlq_500 <- readRDS(here("Objects", "D99_rlq_500.rds"))

#### Custom function - spp ####

# function to create custom RLQ plots for species abundance
# default plot functions in "ade4" have overlapping labels 
# arguments: 
# (1) rlq_obj: "rlq dudi" object
spp_plot <- function(rlq_obj, main_title, y_limit, margins) {
  
   # start of ggplot template
   ggplot() +
  
    # insert vertical line at (0, 0)
    geom_vline(
      xintercept = 0, 
      lwd = 0.5, 
      col = "grey"
      ) + 
     
    # insert horizonal line at (0, 0)
    geom_hline(
      yintercept = 0, 
      lwd = 0.5, 
      col = "grey"
      ) +
    
    # use points to indicate species composition on biplots
    geom_point(
      aes(x = AxcQ1, y = AxcQ2), 
          size = 2, 
          color = "grey", 
          alpha = 1, 
          data = rlq_obj$lQ
      ) +
     
    # create non-overlapping labels from ggrepel package
    # use regular expression to grab genus names for legend
    geom_label_repel(
      aes(x = AxcQ1, 
          y = AxcQ2, 
          fill = (str_extract(rownames(rlq_obj$lQ), "^[^_]+")), #regex, 
          label = rownames(rlq_obj$lQ)),
      segment.alpha = 0.5,
      point.padding = 0.5,
      size = 3,
      data = rlq_obj$lQ
      ) + 
    
    # add in legend title
    scale_fill_discrete(name = "Genus") + 
     
    # increase space for cover type symbols
    ylim(y_limit) + 
     
    # x-axis and y-axis labels
    xlab("Axis 1") + 
    ylab("Axis 2") +
     
    # remove gridlines and grey background
    # add plot margins
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.margin = unit(margins, "cm")
          ) +
    
    # add in titles
    labs(title = main_title)
  
}

# function to create custom RLQ plots for trait states
# default plot functions in "ade4" have overlapping labels 
# arguments: 
# (1) rlq_obj: "rlq dudi" object

#### Custom function - traits ####
trait_plot <- function(rlq_obj, main_title, y_limit, margins) {
  
  #rlq bilots: trait distribution along RLQ axes 
  ggplot() + 
    
  # insert vertical line at (0, 0)
  geom_vline(xintercept = 0, 
             lwd = 0.5, 
             col = "grey"
             ) +  
    
  # insert horizonal line at (0, 0)
  geom_hline(yintercept = 0, 
             lwd = 0.5,
             col = "grey") +
    
  # create non-overlapping labels 
  geom_label_repel(aes(x = CS1, y = CS2, 
                       label = rownames(rlq_obj$c1)), 
                   point.padding = 0.5,
                   data = rlq_obj$c1
                   ) + 
  
  # create arrows from origin (0,0)
  geom_segment(aes(x = 0, y = 0, 
                   xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.01, "npc")),
               colour = "blue",
               alpha = 0.25, 
               data = rlq_obj$c1
               ) +  
                  
  # increase space for cover type symbols
  ylim(y_limit) +
    
  # x-axis and y-axis labels
  xlab("Axis 1") +
  ylab("Axis 2") +
    
  # remove gridlines and grey background 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margins, "cm")
        ) +
  
  # add in titles
  labs(title = main_title)

}

#### Custom function - env ####
env_plot <- function(rlq_obj, main_title, margins) {
    
    # change cover type names
    rownames(rlq_obj$l1) <- c("% Grass", "% Tree Canopy", "% Urban Area", "Edge Density")
    
    # start of ggplot template
    ggplot() +
    
    # insert vertical line at (0, 0)
    geom_vline(xintercept = 0, 
               lwd = 0.5, 
               col = "grey"
    ) +  
    
    # insert horizonal line at (0, 0)
    geom_hline(yintercept = 0, 
               lwd = 0.5,
               col = "grey") +
    
    # create non-overlapping labels 
    geom_label_repel(aes(x = RS1, y = RS2, 
                         label = rownames(rlq_obj$l1)), 
                     point.padding = 0.5,
                     data = rlq_obj$l1 
    ) +
    
    # create arrows from origin (0,0)
    geom_segment(aes(x = 0, y = 0, 
                     xend = RS1, yend = RS2),
                 arrow = arrow(length = unit(0.01, "npc")),
                 colour = "red",
                 size = 1,
                 alpha = 0.25, 
                 data = rlq_obj$l1
    ) +  
    
    # x-axis and y-axis labels
    xlab("Axis 1") +
    ylab("Axis 2") +
    
    # remove gridlines and grey background 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.margin = unit(margins, "cm")
          ) +
    
    # add in titles
    labs(title = main_title)
  
}

#### Custom function - eigs ####
eig_plot <- function(rlq_obj, main_title, margins) {
 
   # create data-frame of eigenvalues and projected inertia
   eigs <- data.frame(Axes = paste("Axis", 1:4, sep = " "),
                      Eigenvalues =  rlq_obj$eig,
                      Projected.Inertia = round(rlq_250$eig/sum(rlq_250$eig)*100, digits = 2)
                     )
  # start plotting
  ggplot(
    data = eigs, 
    aes(x = Axes, y = Projected.Inertia)
    ) +
    
  geom_bar(
    colour = "black", 
    stat = "identity"
      ) + 
    
  # remove gridlines and grey background 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    plot.margin = unit(margins, "cm")
    ) +
    
    # x-axis and y-axis labels
    xlab("Axes") +
    ylab("Projected Inertia (%)") +
    
    # add in titles
    labs(title = main_title)
}


#### Custom function - mult-panel figures ####
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Plots ####

# eigs
eig_250 <- eig_plot(rlq_250, "A)", c(1,0,0,1)) # margins: t, r, b, l
eig_500 <- eig_plot(rlq_500, "A)", c(1,0,0,1)) # margins: t, r, b, l


# env 
env_250 <- env_plot(rlq_250, "B)", c(0, 0, 0, 1)) # margins: t, r, b, l
env_500 <- env_plot(rlq_500, "B)", c(0, 0, 0, 1)) # margins: t, r, b, l


# traits 
trait_250 <- trait_plot(rlq_250, "C)", c(-4, 4), c(0, 1, 1, 0)) # margins: t, r, b, l
trait_500 <- trait_plot(rlq_500, "C)", c(-4, 4), c(0, 1, 1, 0)) # margins: t, r, b, l

# spp
spp_250 <- spp_plot(rlq_250, "D)", c(-5, 5), c(0, 1, 1, 0))
spp_500 <- spp_plot(rlq_500, "D)", c(-5, 5), c(0, 1, 1, 0))


#### Bringing it all together ####

# multi-panel figures
plots_250 <- multiplot(eig_250, spp_250, env_250, trait_250, cols = 2)
plots_500 <- multiplot(eig_500, spp_500, env_500, trait_500, cols = 2)




_

_
