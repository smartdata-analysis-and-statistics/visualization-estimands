library(ggplot2)
library(ggridges) # Joy plots
require(ggVennDiagram)
require(ggpubr)
require(gridExtra)
library(tidyr)
library(dplyr) # Facilitate data manipulation
library(car) # Drawing of ellipses
library(Hmisc)
library(stringr)
library(RColorBrewer)

ven_theme <- function(){
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        plot.title.position = "plot",
        axis.title = element_text(size = 13))
}

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

generate_fig1 <- function() {

  x <- list(C=1:5, B=2:7, A=3:6)
  # "#66C2A5" "#FC8D62" "#8DA0CB" 
  venn <- Venn(x)
  data <- process_data(venn)
  p <- ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
    # 2. set edge layer
    geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 0.7) +
    # 3. set label layer
    # geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
    # 4. region label layer
    scale_color_manual(name = "",
                       values = c("A" = "#66C2A5","B" ="#FC8D62", 'C' = '#8DA0CB'),
                       labels = c('A' = 'Treatment A',
                                  'B' = 'Treatment B',
                                  'C' = 'Treatment C')) +
    xlab(expression("x"[1])) +
    ylab(expression("x"[2])) +
    guides(color = guide_legend(reverse = TRUE))
  
  
  
  p$layers[[1]]$mapping <- aes(fill = name)
  
  yes <- "#E5E7E9"
    no  <- "white"
      isect <- "#D7DBDD"
        
      p1 <- p + scale_fill_manual(values = c(A = yes, 
                                             B..A = isect,
                                             C..A = isect,
                                             C..B..A = isect,
                                             B = yes,
                                             C = yes,
                                             C..B = isect),
                                  guide = "none") +
        ggtitle("A)") +
        ven_theme()
      
      legend <- get_legend(p1)
      
      p1 <- p1 + theme(legend.position="none") 
      
      
      p2 <- p + scale_fill_manual(values = c(A = yes, 
                                             B..A = yes,
                                             C..A = yes,
                                             C..B..A = yes,
                                             B = no,
                                             C = no,
                                             C..B = no),
                                  guide = "none") +
        ggtitle("B)") +
        ven_theme()+ 
        theme(legend.position="none") 
      
      
      p3 <- p + scale_fill_manual(values = c(A = no, 
                                             B..A = yes,
                                             C..A = no,
                                             C..B..A = yes,
                                             B = yes,
                                             C = no,
                                             C..B = yes),
                                  guide = "none") +
        ggtitle("C)") +
        ven_theme()+ 
        theme(legend.position="none") 

      
      p4 <- p + scale_fill_manual(values = c(A = yes, 
                                             B..A = yes,
                                             C..A = no,
                                             C..B..A = yes,
                                             B = no,
                                             C = no,
                                             C..B = no),
                                  guide = "none") +
        ggtitle("D)") +
        ven_theme()+ 
        theme(legend.position="none") 
      
      p5 <- p + scale_fill_manual(values = c(A = yes, 
                                             B..A = no,
                                             C..A = yes,
                                             C..B..A = yes,
                                             B = no,
                                             C = no,
                                             C..B = no),
                                  guide = "none") +
        ggtitle("E)") +
        ven_theme()+ 
        theme(legend.position="none") 
      
      p6 <- p + scale_fill_manual(values = c(A = no, 
                                             B..A = no,
                                             C..A = no,
                                             C..B..A = yes,
                                             B = no,
                                             C = no,
                                             C..B = no),
                                  guide = "none") +
        ggtitle("F)") +
        ven_theme()+ 
        theme(legend.position="none") 

      g <- grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6), ncol = 3)
      ggsave(file="venn.png", width = 20, height = 14, g) #saves g

}

densplot_X <- function(sim, ggtype = "density2d", palette = "Set1") {
  
  set.seed(sim$pars$seed)
  
  dat <- sim$data
  dat$txlabel <- paste("Treatment", dat$treatment)

  g <- ggplot(dat, aes(x = x1, y = x2)) + 
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_line(colour = "grey70", size = 0.2, linetype = 2),
          panel.grid.minor = element_blank()) + 
    xlab("x1") +
    ylab("x2") +
    scale_color_brewer(palette = palette) + 
    scale_fill_brewer(palette = palette)
    
  if (ggtype == "density2d") {
    g <- g + stat_density2d(aes(alpha  = ..level.., color = txlabel), contour = TRUE) +
      guides(alpha = "none")
  } else if (ggtype == "point") {
    means <- dat %>% group_by(txlabel) %>% summarise(mean.x1 = mean(x1), mean.x2 = mean(x2))
    
    if (nrow(dat) > 5000) {
      datsmall <- dat[sample(seq(nrow(dat)), size = 5000),]
    }
    
    g <- g + geom_point(data = datsmall, alpha = 0.4, size = 0.25, aes(color = txlabel)) +
      stat_ellipse(type = "norm", lty = 2, lwd = 1.5,  aes(col = txlabel)) +
      geom_point(means, mapping = aes(x = mean.x1, y = mean.x2, fill = txlabel), size = 4, pch = 21, color = 'black') +
      xlim(-5, 5) +
      ylim(-5, 5)
      
  }
  return(g)
}

plot_joy_x <- function(sim,
                       population,
                       scenario = 1,
                       replace = FALSE,
                       palette = "Set2") {
  
  treatments <- unlist(strsplit(population, split = ""))
  colors <- brewer.pal(3, palette) #A, B, C
  names(colors) <- c("A", "B", "C")
  
  psLabel <- paste("ps", population, sep = "")
  ggdat <- NULL
  if (population == "AB") {
    ggdat <- data.frame(
      "ID" = subset(sim$data, trtC == 0)$ID,
      "treatment" = subset(sim$data, trtC == 0)$treatment,
      "x" = subset(sim$data, trtC == 0)$x1,
      "covariate" = "x1")
    ggdat <- rbind(ggdat, data.frame(
      "ID" = subset(sim$data, trtC == 0)$ID,
      "treatment" = subset(sim$data, trtC == 0)$treatment,
      "x" = subset(sim$data, trtC == 0)$x2,
      "covariate" = "x2"))
  } else if (population == "AC") {
    ggdat <- rbind(ggdat, data.frame(
      "ID" = subset(sim$data, trtB == 0)$ID,
      "treatment" = subset(sim$data, trtB == 0)$treatment,
      "x" = subset(sim$data, trtB == 0)$x1,
      "covariate" = "x1"))
    ggdat <- rbind(ggdat, data.frame(
      "ID" = subset(sim$data, trtB == 0)$ID,
      "treatment" = subset(sim$data, trtB == 0)$treatment,
      "x" = subset(sim$data, trtB == 0)$x2,
      "covariate" = "x2"))
  } else if (population == "BC") {
    ggdat <- rbind(ggdat, data.frame(
      "ID" = subset(sim$data, trtA == 0)$ID,
      "treatment" = subset(sim$data, trtA == 0)$treatment,
      "x" = subset(sim$data, trtA == 0)$x1,
      "covariate" = "x1"))
    ggdat <- rbind(ggdat, data.frame(
      "ID" = subset(sim$data, trtA == 0)$ID,
      "treatment" = subset(sim$data, trtA == 0)$treatment,
      "x" = subset(sim$data, trtA == 0)$x2,
      "covariate" = "x2"))
  } else {
    stop("Invalid population!")
 }
  
  ggdat$wtype <- "Original sample"
  
  
  for (fname in retrieve_ATT_files(population, replace, fdir = paste0("Datasets/Scenario ", scenario, "/"))){
    fit <- NULL
    load(paste(paste0("Datasets/Scenario ", scenario, "/"), fname, sep = ""))
    dati <- get_matches(fit)[,c("ID", "treatment", "x1", "x2")]
    dati <- data.frame(dati, wtype = paste("ATT-", fit$target_population, sep = ""))
    ggdat <- rbind(ggdat, data.frame(
      "ID" = dati$ID,
      "treatment"= dati$treatment,
      "x" = dati$x1,
      "covariate" = "x1",
      wtype = paste("ATT-", fit$target_population, sep = "")
    ))
    ggdat <- rbind(ggdat, data.frame(
      "ID" = dati$ID,
      "treatment"= dati$treatment,
      "x" = dati$x2,
      "covariate" = "x2",
      wtype = paste("ATT-", fit$target_population, sep = "")
    ))
  }

  ggdat$wtype <- factor(ggdat$wtype, levels = rev(unique(ggdat$wtype)))
  
  vline.data <- ggdat %>%
                    group_by(wtype, covariate,treatment) %>%
                    dplyr::summarize(z = mean(x))
  
  ggplot(ggdat) + 
    geom_density(aes(x = x, fill = treatment, color = treatment),
                 alpha = 0.25)+

    scale_color_manual(values = colors[treatments],
                       labels = c(paste("Treatment", treatments[1]), paste("Treatment", treatments[2]))) + 
    scale_fill_manual(values = colors[treatments],
                      labels = c(paste("Treatment", treatments[1]), paste("Treatment", treatments[2])))+
    #labs(title = "Distribution of the Propensity Score") +
    ylab("Density")+ 
    xlab("")+
    theme(strip.placement = "outside") +
    facet_grid(wtype~covariate, 
               scales = "free_x", switch = 'x') +
    theme(axis.title.x = element_blank(),
          strip.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank())+
    geom_vline(aes(xintercept = z, color = treatment, linetype = treatment), vline.data, size = 1.2)+
    scale_linetype_manual(values=c("dashed", "dotted"), guide = "none")
}

plot_distr_ps <- function(sim, 
                          population, # Which population to extract? "AB", "BC" or "AC"
                          replace = FALSE,
                          palette = "Set1") 
  {
  treatments <- unlist(strsplit(population, split = ""))
  colors <- brewer.pal(3, palette) #A, B, C
  names(colors) <- c("A", "B", "C")
  
  psLabel <- paste("ps", population, sep = "")
  
  if (population == "AB") {
    ggdat <- subset(sim$data, trtC == 0)[,c("ID", "treatment",psLabel)]
  } else if (population == "AC") {
    ggdat <- subset(sim$data, trtB == 0)[,c("ID", "treatment",psLabel)]
  } else if (population == "BC") {
    ggdat <- subset(sim$data, trtA == 0)[,c("ID", "treatment",psLabel)]
  } else {
    stop("Invalid population!")
  }
  ggdat$wtype <- "Original sample"
  
  for (fname in retrieve_ATT_files(population, replace, fdir = sim$pars$save.dir)) {
    fit <- NULL
    load(paste(sim$pars$save.dir, fname, sep = ""))
    dati <- get_matches(fit)[,c("ID", "treatment", psLabel)]
    dati <- data.frame(dati, wtype = paste("ATT", fit$target_population, sep = ""))
    ggdat <- rbind(ggdat, dati)
  }
  ggdat$uniqrow <- seq(nrow(ggdat))
  ggdat$wtype <- factor(ggdat$wtype, levels = unique(ggdat$wtype), labels = unique(ggdat$wtype))
  
  d <- ggdat %>%
    spread(treatment, paste(psLabel), sep = "_p")
  
  lbl_treatment1 <- paste("treatment_p", treatments[1], sep = "")
  lbl_treatment2 <- paste("treatment_p", treatments[2], sep = "")
  
  
  
  
  ggplot(d) + 
    geom_histogram(breaks = seq(0, 1, length = 50), aes(x = d[,lbl_treatment1]),
                   fill = colors[treatments[1]],
                   color = colors[treatments[1]],
                   alpha = 0.5) + 
    geom_histogram(breaks = seq(0, 1, length = 50), aes(x = d[,lbl_treatment2], y = -..count..),
                   fill = colors[treatments[2]],
                   color = colors[treatments[2]],
                   alpha = 0.5) +
    ylab("Number of patients") + 
    xlab("Probability of receiving drug A") +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = abs) +
    labs(title = "Distribution of the Propensity Score") +
    theme(legend.position = "bottom") +
    xlim(c(0,1)) +
    facet_wrap(~wtype) +
    scale_color_manual(name = 'the colour', 
                       values = colors[treatments],
                       labels = c(paste("Patients receiving", treatments[1]), paste("Patients receiving", treatments[2]))) + 
    scale_fill_manual(values = colors[treatments],
                      labels = c(paste("Patients receiving", treatments[1]), paste("Patients receiving", treatments[2])))
}

prepare_ggdat_original_sample <- function(sim) {
  
  ggdat <- NULL
  
  ##############################################################################
  # Prepare AB population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtC == 0) %>% select(ID),
    "x1" = (sim$data %>% filter(trtC == 0))$x1,
    "x2" = (sim$data %>% filter(trtC == 0))$x2,
    "treatment" = (sim$data %>% filter(trtC == 0))$treatment,
    "ps" = (sim$data %>% filter(trtC == 0))$psAB,
    "population" = "AB population",
    "wtype" = "Original sample"
  ))
  
  ##############################################################################
  # Prepare AC population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtB == 0) %>% select(ID),
    "x1" = (sim$data %>% filter(trtB == 0))$x1,
    "x2" = (sim$data %>% filter(trtB == 0))$x2,
    "treatment" = (sim$data %>% filter(trtB == 0))$treatment,
    "ps" = (sim$data %>% filter(trtB == 0))$psAC,
    "population" = "AC population",
    "wtype" = "Original sample"
  ))
  
  ##############################################################################
  # Prepare BC population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtA == 0) %>% select(ID),
    "x1" = (sim$data %>% filter(trtA == 0))$x1,
    "x2" = (sim$data %>% filter(trtA == 0))$x2,
    "treatment" = (sim$data %>% filter(trtA == 0))$treatment,
    "ps" = (sim$data %>% filter(trtA == 0))$psBC,
    "population" = "BC population",
    "wtype" = "Original sample"
  ))
  
  ggdat
}

plot_all_distr_ps <- function(sim, 
                          replace = FALSE,
                          palette = "Set2") # show results for matching with replacement?
{
  treatments <- c("A", "B", "C")
  colors <- brewer.pal(3, palette) #A, B, C
  names(colors) <- c("A", "B", "C")
  
 
  ggdat <- NULL
  
  ##############################################################################
  # Prepare AB population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtC == 0) %>% select(ID),
    "treatment" = (sim$data %>% filter(trtC == 0))$treatment,
    "ps" = (sim$data %>% filter(trtC == 0))$psAB,
    "population" = "AB population",
    "wtype" = "Original sample"
  ))
  
  ##############################################################################
  # Prepare AC population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtB == 0) %>% select(ID),
    "treatment" = (sim$data %>% filter(trtB == 0))$treatment,
    "ps" = (sim$data %>% filter(trtB == 0))$psAC,
    "population" = "AC population",
    "wtype" = "Original sample"
  ))
  
  ##############################################################################
  # Prepare BC population
  ##############################################################################
  ggdat <- rbind(ggdat, data.frame(
    "ID" = sim$data %>% filter(trtA == 0) %>% select(ID),
    "treatment" = (sim$data %>% filter(trtA == 0))$treatment,
    "ps" = (sim$data %>% filter(trtA == 0))$psBC,
    "population" = "BC population",
    "wtype" = "Original sample"
  ))
    
  for (population in c("AB", "AC", "BC")) {
    for (fname in retrieve_ATT_files(population, replace, fdir ="Datasets/Scenario 1/")) {
      fit <- NULL
      load(paste("Datasets/Scenario 1/", fname, sep = ""))
      ggdat <- rbind(ggdat, data.frame(
        "ID" = get_matches(fit)$ID,
        "treatment" = get_matches(fit)$treatment,
        "ps" = get_matches(fit)[,paste("ps", population, sep = "")],
        "population" = paste(population, "population"),
        "wtype" = paste("ATT", fit$target_population, sep = "")
      ))
    }
  }
  
  
  ggdat$uniqrow <- seq(nrow(ggdat))
  ggdat$wtype <- factor(ggdat$wtype, levels = unique(ggdat$wtype), labels = unique(ggdat$wtype))
  ggdat$treatment_p1 <- ggdat$treatment_p2 <- NA
  
  ## AB comparisons
  ind_ATTA <- which(ggdat$population == "AB population" & ggdat$treatment == "A")
  ind_ATTB <- which(ggdat$population == "AB population" & ggdat$treatment == "B")
  ggdat$treatment_p1[ind_ATTA] <- ggdat$ps[ind_ATTA]
  ggdat$treatment_p2[ind_ATTB] <- ggdat$ps[ind_ATTB]
  
  ## AC comparisons
  ind_ATTA <- which(ggdat$population == "AC population" & ggdat$treatment == "A")
  ind_ATTC <- which(ggdat$population == "AC population" & ggdat$treatment == "C")
  ggdat$treatment_p1[ind_ATTA] <- ggdat$ps[ind_ATTA]
  ggdat$treatment_p2[ind_ATTC] <- ggdat$ps[ind_ATTC]

  ## BC comparisons
  ind_ATTB <- which(ggdat$population == "BC population" & ggdat$treatment == "B")
  ind_ATTC <- which(ggdat$population == "BC population" & ggdat$treatment == "C")
  ggdat$treatment_p1[ind_ATTB] <- ggdat$ps[ind_ATTB]
  ggdat$treatment_p2[ind_ATTC] <- ggdat$ps[ind_ATTC]
  
  
  # draw grey rectangle
  rect <- data.frame(population = c("BC population", "AC population", "AB population"),
                     wtype = c("ATTA", "ATTB", "ATTC"),
                     xmin = c(0, 0, 0), 
                     xmax = c(1, 1, 1), 
                     ymin = c(-3000, -3000, -3000), 
                     ymax = c(3000, 3000, 3000))
  rect$population <- factor(rect$population, levels = unique(ggdat$population), labels = unique(ggdat$population))
  rect$wtype <- factor(rect$wtype, levels = unique(ggdat$wtype), labels = unique(ggdat$wtype))
  
  
  
  ggplot(ggdat) + 
    geom_histogram(breaks = seq(0, 1, length = 50), 
                   aes(x = treatment_p1, fill = treatment, color = treatment),
                   alpha = 0.5) + 
    geom_histogram(breaks = seq(0, 1, length = 50), 
                   aes(x = treatment_p2, y = -..count.., fill = treatment, color = treatment),
                   alpha = 0.5) +
    ylab("Number of patients") + 
    xlab("Propensity Score") +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = abs) +
    labs(title = "Distribution of the Propensity Score") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_line(colour = "grey70", size = 0.2, linetype = 2),
          panel.grid.minor = element_blank()) +
    xlim(c(0,1)) +
    facet_grid(population ~ wtype) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              alpha = 1, fill = "#D3D3D3",
              data = rect,
              inherit.aes = FALSE)+
    scale_color_manual(name = 'the colour', 
                       values = colors,
                       labels = paste("Patients receiving", treatments),
                       guide = "none") + 
    scale_fill_manual(name = "", 
                      values = colors,
                      labels = paste("Patients receiving", treatments))
}

densplot_all_matched_samples <- function(sim, 
                                         replace = FALSE, 
                                         ggtype = "density2d",
                                         palette = "Set1") 
{
  treatments <- c("A", "B", "C")
  colors <- brewer.pal(3, palette) #A, B, C
  txlabels <- names(colors) <- c("Treatment A", "Treatment B", "Treatment C")
  fdir <- "Datasets/Scenario 1/"
  
  ggdat <- prepare_ggdat_original_sample(sim)
  
  for (population in c("AB", "AC", "BC")) {
    for (fname in retrieve_ATT_files(population, replace, fdir = fdir)) {
      fit <- NULL
      load(paste0(fdir, fname))
      ggdat <- rbind(ggdat, data.frame(
        "ID" = get_matches(fit)$ID,
        "x1" = get_matches(fit)$x1,
        "x2" = get_matches(fit)$x2,
        "treatment" = get_matches(fit)$treatment,
        "ps" = get_matches(fit)[,paste0("ps", population)],
        "population" = paste(population, "population"),
        "wtype" = paste0("ATT", fit$target_population)
      ))
    }
  }
  
  ggdat$population <- factor(ggdat$population, levels = c("AB population", "AC population", "BC population"), 
                             labels = c("A-B comparison", "A-C comparison", "B-C comparison"))
  ggdat$wtype <- factor(ggdat$wtype, levels =  c("Original sample", "ATTA", "ATTB", "ATTC"), 
                        labels = c("Original sample", "Matched ATT-A", "Matched ATT-B", "Matched ATT-C"))
  ggdat$txlabel <- paste("Treatment", ggdat$treatment)
  
  rect <- data.frame(population = c("B-C comparison", "A-C comparison", "A-B comparison", "A-B comparison"),
                   wtype = c("Matched ATT-A", "Matched ATT-B", "Matched ATT-C", "Original sample"),
                   xmin = c(-5, -5, -5, -10), 
                   xmax = c(5, 5, 5, -9), 
                   ymin = c(-5, -5, -5, -10), 
                   ymax = c(5, 5, 5, -9))
  rect$population <- factor(rect$population, levels = unique(ggdat$population), labels = unique(ggdat$population))
  rect$wtype <- factor(rect$wtype, levels = unique(ggdat$wtype), labels = unique(ggdat$wtype))
  
  
  g <- ggplot(ggdat, aes(x = x1, y = x2)) +
    facet_grid(population ~ wtype) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              alpha = 1, fill = "#D3D3D3",
              data = rect,
              inherit.aes = FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_line(colour = "grey70", size = 0.2, linetype = 2),
          panel.grid.minor = element_blank()) +
    scale_color_manual(values = colors, labels = txlabels) + 
    scale_fill_manual(values = colors, labels = txlabels)
  
  if (ggtype == "density2d") {
    g <- g + stat_density2d(aes(alpha  = ..level.., color = txlabel), contour = TRUE) +
      guides(alpha = "none") + # Dont include alpha level in legend 
      labs(title = "Covariate Overlap Pre- and Post-Matching",
           subtitle = "Bivariate Density Contours")
  } else if (ggtype == "point") {
    # get ellipse data with and without weighting
    edata1 <- ggdat %>%
      group_by(txlabel, wtype, population) %>%
      do({
        
        x <- .$x1
        y <- .$x2
        # get boundary data for bivariate 95% confidence ellipse 
        out <- as.data.frame(dataEllipse(x, y, draw = FALSE)$`0.95`)
        # rename columns to original 
        colnames(out) <- c("x1", "x2")
        
        out
      })
    
    sdata <- ggdat %>%
      group_by(txlabel, wtype, population) %>%
      summarise(mean_x1 = mean(x1),
                mean_x2 = mean(x2))
    
    if (nrow(ggdat) > 5000) {
      ggdatsmall <- ggdat[sample(seq(nrow(ggdat)), size = 5000),]
    }
    
    g <- g + geom_point(data = ggdatsmall, 
                        alpha = 0.4, 
                        aes(color = txlabel), 
                        size = 0.25) +
      geom_polygon(data = edata1, 
                   fill = NA,
                   size = 1, 
                   aes(linetype = txlabel, color = txlabel)) +
      scale_linetype_manual(values = c("Treatment A" = "dotted", 
                                       "Treatment B" = "dashed",
                                       "Treatment C" = "longdash")) +
      geom_point(data = sdata, aes(mean_x1, mean_x2, fill = txlabel), 
                 size = 4, 
                 alpha = 0.5, pch = 21, color = 'black') +
      ylim(-5,5) +
      xlim(-5,5)
  }
  
  return(g)
}


densplot_matched_samples <- function(sim, 
                                     population, # Which population to extract? "AB", "BC" or "AC"
                                     replace = FALSE, # show results for matching with replacement?
                                     ggtype = "density2d",
                                     palette = "Set1") 
{
  treatments <- unlist(strsplit(population, split = ""))
  colors <- brewer.pal(3, palette) #A, B, C
  names(colors) <- c("A", "B", "C")
  
  if (population == "AB") {
    data <- subset(sim$data, trtC == 0)
  } else if (population == "AC") {
    data <- subset(sim$data, trtB == 0)
  } else if (population == "BC") {
    data <- subset(sim$data, trtA == 0)
  } else {
    stop("Invalid population!")
  }
  
  
  
  ggdat <- data.frame("ID" = character(),
                      "x1" = numeric(),
                      "x2" = numeric(),
                      "treatment" = character(),
                      "wtype" = character())
  
  ## Add info from the original sample
  ggdat <- rbind(ggdat, data.frame("ID" = data$ID,
                                   "x1" = data$x1,
                                   "x2" = data$x2,
                                   "treatment" = as.character(data$treatment),
                                   "wtype" = "Original sample"))
    
  for (fname in retrieve_ATT_files(population, replace, fdir = sim$pars$save.dir)) {
    fit <- NULL
    load(paste(sim$pars$save.dir, fname, sep = ""))
    mdati <-  get_matches(fit, data = data)
    ggdat <- rbind(ggdat, data.frame("ID" = mdati$ID,
                                     "x1" = mdati$x1,
                                     "x2" = mdati$x2,
                                     "treatment" = as.character(mdati$treatment),
                                     "wtype" = paste("ATT", fit$target_population, sep = "")))
  }

  ggdat$wtype <- factor(ggdat$wtype, levels = unique(ggdat$wtype), labels = unique(ggdat$wtype))
  
 
  
  g <- ggplot(ggdat, aes(x = x1, y = x2)) +
    facet_wrap(~wtype) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = colors[treatments], labels = treatments) + 
    scale_fill_manual(values = colors[treatments], labels = treatments)
  
  if (ggtype == "density2d") {
    g <- g + stat_density2d(aes(alpha  = ..level.., color = treatment), contour = TRUE) +
      guides(alpha = "none") + # Dont include alpha level in legend 
      labs(title = "Covariate Overlap Pre- and Post-Matching",
           subtitle = "Bivariate Density Contours")
  } else if (ggtype == "point") {
    # get ellipse data with and without weighting
    edata1 <- ggdat %>%
      group_by(treatment, wtype) %>%
      do({
        
        x <- .$x1
        y <- .$x2
        
        # get boundary data for bivariate 95% confidence ellipse 
        out <- as.data.frame(dataEllipse(x, y, draw = FALSE)$`0.95`)
        # rename columns to original 
        colnames(out) <- c("x1", "x2")
        
        out
      })
    
    sdata <- ggdat %>%
      group_by(treatment, wtype) %>%
      summarise(mean_x1 = mean(x1),
                mean_x2 = mean(x2))
    
    if (nrow(ggdat) > 5000) {
      ggdatsmall <- ggdat[sample(seq(nrow(ggdat)), size = 5000),]
    }
    
    g <- g + geom_point(data = ggdatsmall, alpha = 0.4, aes(color = treatment), size = 0.25) +
      geom_polygon(data = edata1, fill = NA, size = 2, linetype = "dashed", aes(color = treatment)) +
      geom_point(data = sdata, aes(mean_x1, mean_x2, fill = treatment), size = 4, alpha = 1, pch = 21, color = 'black') +
      labs(title = "Covariate Overlap Pre- and Post-Matching",
           subtitle = "Bivariate 95% Confidence Regions")
  }
  
  return(g)
}

save_all_loveplots <- function(scenario, saveimg = T){
  library(gridExtra)
  library(ggpubr)
  library(cobalt)  
  
  files <- c("A_AB",
             "A_AC",
             "B_AB",
             "B_BC",
             "C_AC",
             "C_BC")
  plot_loveplot <- function(fit){
    plot <- love.plot(fit, 
                      drop.distance = TRUE, 
                      var.order = "unadjusted",
                      abs = TRUE,
                      line = TRUE, stars = "raw",
                      sample.names = c("Pre-matching", "Post-matching"),
                      title = "", 
                      thresholds = c(m = .1)) + theme(legend.title=element_blank())
    return(plot)
  }
  plots <- list()
  for(i in files){
    #cat("Plotting", i, "\n")
    fit <- NULL
    load(paste0("Datasets/", scenario, "/ATT", i, "_with_replacement.rda"))
    plots[[i]] <- plot_loveplot(fit)
  }
  
  if(saveimg == T) {
  ggarrange(plotlist = plots, nrow = 3, ncol = 2,  labels = c("A)", "B)", "C)", "D)", "E)", "F)"), common.legend = TRUE, legend="bottom") %>% 
    ggexport(filename = paste0("Figures/", scenario, "/All love plots.png"),width = 900, height = 900, res = 100)
  } else {
    ggarrange(plotlist = plots, nrow = 3, ncol = 2,  labels = c("A)", "B)", "C)", "D)", "E)", "F)"), common.legend = TRUE, legend="bottom")
  }
}


save_all_balanceplots <- function(scenario, saveimg = T){
  plots <- list()
  colors1 <- c(2,1,3,1,3,2)
  colors2 <- c(1,2,1,3,2,3)
  files <- c("A_AB",
             "A_AC",
             "B_AB",
             "B_BC",
             "C_AC",
             "C_BC")
  colors <- brewer.pal(3, "Set2")
  count <- 0
  for(i in files){
    #cat("Plotting", i, "\n")
    count <- count + 1
    fit <- NULL
    load(paste0("Datasets/", scenario,"/ATT", i, "_with_replacement.rda"))
    plots[[i]] <- bal.plot(fit, var.name = "distance", which = "both",
                           type = "histogram", mirror = T, sample.names = c("Pre-matching", "Post-matching")) + 
      scale_fill_manual(name = "", values = c(colors[colors1[count]], colors[colors2[count]]), 
                        labels = c(ifelse(fit$target_population == str_split(fit$target_comparison, "")[[1]][2], str_split(fit$target_comparison, "")[[1]][1], 
                                          str_split(fit$target_comparison,"")[[1]][2]),  fit$target_population))+
      scale_x_continuous(name= paste("Probability of being treated with", fit$target_population))+
      scale_y_continuous(limits = c(-0.27, 0.27), breaks = c(-0.2, -0.1, 0, 0.1, 0.2), labels= c(0.2, 0.1, 0, 0.1, 0.2))+
      ggtitle("")
    
  }
  
  if(saveimg == T) {
    ggarrange(plotlist = plots, nrow = 3, ncol = 2,  labels = c("A)", "B)", "C)", "D)", "E)", "F)"), common.legend = F, legend="bottom") %>%
              ggexport(filename = paste0("Figures/", scenario, "/All balance plots.png"), width = 900, height = 900, res = 100) 
  } else {
    ggarrange(plotlist = plots, nrow = 3, ncol = 2,  labels = c("A)", "B)", "C)", "D)", "E)", "F)"), common.legend = F, legend="bottom") 
  }
}

generate_figures <- function(scenario_num = 1) {
  scenario_lbl <- paste("Scenario", scenario_num)
  suppressMessages({
  # Load simulation study results 
  load(paste0("Datasets/", scenario_lbl, "/simdat.rda")) # Simulated data
  load(paste0("Datasets/", scenario_lbl, "/results.rda")) # ATT estimates
  
  ################################################################################
  # Ellipses graph
  ################################################################################
  densplot_X(sim, ggtype = "point", palette = "Set2")
  ggsave("Figures/Scenario 1/ellipses.png")
  
  ################################################################################
  # Distribution of the propensity score [Panel of all 3 plots]
  ################################################################################
  plot_distr_ps(sim, population = "AB", replace = TRUE)
  ggsave(paste0("Figures/", scenario_lbl, "/ps_distr_AB_matching_with_replacement.png"))
  plot_distr_ps(sim, population = "BC", replace = TRUE)
  ggsave(paste0("Figures/", scenario_lbl, "/ps_distr_BC_matching_with_replacement.png"))
  plot_distr_ps(sim, population = "AC", replace = TRUE)
  ggsave(paste0("Figures/", scenario_lbl, "/ps_distr_AC_matching_with_replacement.png"))
  
  ################################################################################
  # Covariate overlap [Panel of all 3 plots]
  ################################################################################
  densplot_matched_samples(sim, population = "AB", replace = TRUE, ggtype = "point")
  ggsave(paste0("Figures/", scenario_lbl, "/covl_AB_matching_with_replacement.png"))
  densplot_matched_samples(sim, population = "BC", replace = TRUE, ggtype = "point")
  ggsave(paste0("Figures/", scenario_lbl, "/covl_BC_matching_with_replacement.png"))
  densplot_matched_samples(sim, population = "AC", replace = TRUE, ggtype = "point")
  ggsave(paste0("Figures/", scenario_lbl, "/covl_AC_matching_with_replacement.png"))
  
  ################################################################################
  # Joy plot for only ATTA and ATTB
  ################################################################################
  plot_joy_x(sim, population = "AB", replace = TRUE)
  ggsave(paste0("Figures/", scenario_lbl, "/joy_AB_matching_with_replacement.png"))
  
  ################################################################################
  # Panel plot of all "Distribution of the propensity scores" graphs
  ################################################################################
  plot_all_distr_ps(sim, replace = TRUE)
  ggsave(paste0("Figures/", scenario_lbl, "/ps_distr_all_matching_with_replacement.png"), width = 15, height = 12, scale = 0.6) # SW: added scale to increase font size
  
  ################################################################################
  # Panel plot of all"Covariate overlap Pre- and Post-matching (bivariate 95% Confidence Regions only)" graphs
  ################################################################################
  densplot_all_matched_samples(sim, replace = TRUE, ggtype = "point", palette = "Set2")
  ggsave(paste0("Figures/", scenario_lbl, "/covl_all_matching_with_replacement.png"))
  
  ################################################################################
  # Combined Love plots and balance plots
  ################################################################################
  save_all_loveplots(scenario = scenario_lbl)
  save_all_balanceplots(scenario = scenario_lbl)
  })
   
  cat(paste0("All plots for ", scenario_lbl, " are saved in the /Figures/", scenario_lbl, " folder."))
}

