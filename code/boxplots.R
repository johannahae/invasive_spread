
# Functions to generate boxplots

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)) #95%CI
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

abiotic.boxplot <- function(dat=datin, x=x, xu=xu, y=y, pal=pal, 
                            namey="Fraction successful invasions"){
  p1 <- dat %>%  
  ggplot(aes(x=as.factor(x), y=y, alpha=0.4)) + 
  stat_summary(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4, 
               show.legend = FALSE, position = position_dodge(1),
               # show.legend = TRUE, position = position_dodge(1),
               fun.data = quantiles_95, geom="boxplot") +
  scale_y_continuous(name =  paste0(namey, "\n"),
                     breaks = seq(0, 100, 50),
                     limits=c(0, 105), 
                     labels = function(b) { paste0(round(b * 1, 0), "%")}) +  
  scale_x_discrete(name = "",
                   breaks = xu,
                   limits = xu, 
                   labels = xu) +
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap(~nutrient.scenario, labeller=label_value) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size=28), 
        axis.text.x=element_blank(),
        axis.title.y.left = element_text(size = 44), 
        axis.text.y.left = element_text(size = 28),  
        axis.ticks.y.left = element_line(size=1),  
        axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 40), 
        legend.title = element_blank(), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width  = unit(1, "cm"), 
        legend.text = element_text(size=36), 
        legend.spacing.y = unit(0.25, "cm"), 
        legend.position = "bottom") +
  guides(alpha=FALSE, fill = guide_legend(reverse=FALSE, override.aes = list(alpha = 0.6)),
         col = guide_legend(reverse=FALSE))
  
  p1 + stat_compare_means(comparisons = list(c("random", "clustered")), 
                          label = "p.signif", method = "t.test")
}

abiotic.boxplot.allnuts <- function(dat=datin, x=x, xu=xu, y=y, pal=pal, 
                                     namey="Fraction successful invasions"){
  p1 <- dat %>%  
    ggplot(aes(x=as.factor(x), y=y, alpha=0.4)) + 
    stat_summary(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4,
                 show.legend = FALSE, position = position_dodge(1),
                 fun.data = quantiles_95, geom="boxplot") +
    scale_y_continuous(name =  paste0(namey, "\n"),
                       breaks = seq(0, 100, 50),
                       limits=c(0, 105), 
                       labels = function(b) { paste0(round(b * 1, 0), "%")}) +  
    scale_x_discrete(name = "",
                     breaks = xu,
                     limits = xu, 
                     labels = xu) +
    scale_colour_manual(values=pal) +
    scale_fill_manual(values=pal) +
    ggpubr::theme_pubr() +
    theme(text = element_text(size=28), 
          axis.text.x=element_blank(),
          axis.title.y.left = element_text(size = 44), 
          axis.text.y.left = element_text(size = 28),  
          axis.ticks.y.left = element_line(size=1),  
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 40), 
          legend.title = element_blank(), 
          legend.key.height = unit(2, "cm"), 
          legend.key.width  = unit(1, "cm"), 
          legend.text = element_text(size=36), 
          legend.spacing.y = unit(0.25, "cm"), 
          legend.position = "bottom") +
    guides(alpha=FALSE, fill = guide_legend(reverse=FALSE, override.aes = list(alpha = 0.6)),
           col = guide_legend(reverse=FALSE))
    
    p1 + stat_compare_means(comparisons = list(c("random", "clustered")), 
                            label = "p.signif", method = "t.test")
}

abiotic.boxplot.log <- function(dat=datin, x=x, xu=xu, y=y, pal=pal, 
                                namey="Fraction successful invasions"){
  p1 <- dat %>% 
    ggplot(aes(x=as.factor(x), y=y, alpha=0.5)) + 
    stat_summary(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4, 
                 show.legend = FALSE, position = position_nudge(x=-0.0), 
                 fun.data = quantiles_95, geom="boxplot") +
    scale_y_log10(name = paste0(namey, "\n"), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(name = "") +
    scale_colour_manual(values=pal) +
    scale_fill_manual(values=pal) +
    facet_wrap(~nutrient.scenario, labeller=labeller(label_value)) + 
    ggpubr::theme_pubr() +
    theme(text = element_text(size=28), 
          # axis.text.x=element_text(angle = 45, hjust = 1 , colour = pal),
          axis.text.x=element_blank(),         
          axis.title.y.left = element_text(size = 44), 
          axis.text.y.left = element_text(size = 28),  
          axis.ticks.y.left = element_line(size=1),  
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 40), 
          legend.title = element_blank(), 
          legend.key.height = unit(2, "cm"), 
          legend.key.width  = unit(1, "cm"), 
          legend.text = element_text(size=36), 
          legend.spacing.y = unit(0.25, "cm"), 
          legend.position = "bottom") +
    guides(alpha=FALSE, fill = guide_legend(reverse=FALSE, override.aes = list(alpha = 0.6)),
           col = guide_legend(reverse=FALSE))
  
  a <- annotation_logticks(base=10, sides="l", scaled = TRUE)
  a$data <- data.frame(x=NA, nutrient.scenario="oligotrophic")
  
  p1 + a + stat_compare_means(comparisons = list(c("random", "clustered")), 
                          label = "p.signif", method = "t.test")
}


biotic.boxplot.right <- function(dat=datin, x=x, xu=xu, pal=pal, barwidth = w, 
                            namex="Dispersal distance", xmin=-0.5, xmax=5, mylabels=xu/10) {
  dat %>% ggplot(aes(x=x, y=fraction.of.invaded.patches), alpha=1) + 
    geom_boxplot(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4, width = barwidth,
                 show.legend = FALSE, position = position_dodge2(preserve = "single"), notch=FALSE,
                 stat = "summary", fun.data = quantiles_95) +
    scale_y_continuous(name = "Fraction invaded patches\n",
                       breaks = seq(0, 100, 50),
                       limits=c(0, 100), 
                       labels = function(b) { paste0(round(b * 1, 0), "%")}) + 
    scale_x_discrete(name = paste("\n", namex),
                     breaks = xu,
                     limits = xu, 
                     labels = mylabels) +
    scale_colour_manual(values=pal) +
    scale_fill_manual(values=pal) + 
    facet_grid(landscape.type ~ nutrient.scenario, 
               labeller=labeller(.rows = label_value, .cols = label_value)) + 
    ggpubr::theme_pubr()+
    expand_limits(x=c(xmin,xmax)) +
    theme(text = element_text(size=20), 
          strip.text.y = element_text(size=40),
          # strip.text.y = element_blank(),
          strip.text.x = element_text(size=40),
          axis.ticks.x = element_line(size = 1),
          axis.ticks.y = element_line(size = 1),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 44))
}

biotic.boxplot.left <- function(dat=datin, x=x, xu=xu, pal=pal, barwidth = w, 
                           namex="Dispersal distance", xmin=-0.5, xmax=5, mylabels=xu/10) {
  dat %>% ggplot(aes(x=x, y=fraction.of.invaded.patches), alpha=1) + 
    geom_boxplot(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4, width = barwidth,
                 show.legend = FALSE, position = position_dodge2(preserve = "single"), notch=FALSE,
                 stat = "summary", fun.data = quantiles_95) +
    scale_y_continuous(name = "Fraction invaded patches\n",
                       breaks = seq(0, 100, 50),
                       limits=c(0, 100), 
                       labels = function(b) { paste0(round(b * 1, 0), "%")}) + 
    scale_x_discrete(name = paste("\n", namex),
                     breaks = xu,
                     limits = xu, 
                     labels = mylabels) +
    scale_colour_manual(values=pal) +
    scale_fill_manual(values=pal) + 
    facet_grid(landscape.type ~ nutrient.scenario, 
               labeller=labeller(.rows = label_value, .cols = label_value)) + 
    ggpubr::theme_pubr() +
    expand_limits(x=c(xmin,xmax)) +
    theme(text = element_text(size=20), 
          strip.text.y = element_blank(),
          strip.text.x = element_text(size=40),
          axis.ticks.x = element_line(size = 1),
          axis.ticks.y = element_line(size = 1),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          axis.title.y = element_text(size = 44),
          axis.title.x = element_text(size = 44))
}