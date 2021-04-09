

#### Declare libraries ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(igraph)
library(NetIndices)
library(RColorBrewer)
library(here)

# Declare path
here::i_am("invasive_spread.Rproj")
# here::here()

# Define and set color palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7") # palette with grey
gg_color_hue <- function(n) cbPalette[1:n]

pal <- gg_color_hue(8)
ramp <- colorRampPalette(pal,bias=0.75)
spcols <- c(ramp(21)[3], ramp(21)[19])

bg <- c("blue", "darkgreen")

# Define and set theme
theme_half_open(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 11/14, rel_tiny = 12/14, rel_large = 16/14)

theme_set(theme_half_open() + background_grid())

# Import data and functions
dat_plot <- here("data/results/dat_plot.rds") %>% read_rds()
source(here("code/boxplots.R"))


#### Figure S1A: Trophic level t = 0 ~ log10 body mass ###
figS1a <- ggscatter(dat_plot, x = "log10mass", y = "trophic.level.invspp",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Body mass [log10]", ylab ="Trophic level at t=0",
                color="invspp.type", palette=spcols, facet.by = "invspp.type",
                panel.labs = list(invspp.type=c("Plant species", "Animal species")),
                legend.title=" ") 
# figS1a %>% ggsave(filename=here("figures/corr_log10mass_trophiclevel0.pdf"),width=16,height=8)

#### Figure S1B: Trophic level t = 5.0000 ~ log10 body mass ###
figS1b <- ggscatter(dat_plot, x = "log10mass", y = "trophic.level.invspp.int",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Body mass [log10]", ylab = "Trophic level at t=5.000",
                color="invspp.type", palette=spcols, facet.by = "invspp.type",
                panel.labs = list(invspp.type=c("Plant species", "Animal species")), 
                legend.title=" ") 
# figS1b %>% ggsave(filename=here("figures/corr_log10mass_trophiclevel5000.pdf"),width=16,height=8)

#### Figure S2: Max. dispersal distance ~ log10 body mass ###
figS2 <- ggscatter(dat_plot, x = "log10mass", y = "dispersal.dist.invspp",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Body mass [log10]", ylab = "Maximum dispersal distance", 
          color="invspp.type", palette=spcols, facet.by = "invspp.type",
          panel.labs = list(invspp.type=c("Plant species", "Animal species")),
          legend.title=" ") 
# figS2 %>% ggsave(filename=here("figures/corr_log10mass_maxdispersalrange.pdf"),width=16,height=8)

#### Figure S3: Realized landscape connectance ~ max. dispersal distance ###
figS3 <- ggscatter(dat_plot, x = "dispersal.dist.invspp", y = "rgg.con.invspp",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Maixmum dispersal distance", ylab = "Realized landscape connectance", 
          color="invspp.type", palette=spcols, 
          facet.by = c("landscape.type","invspp.type"),
          panel.labs = list(invspp.type=c("Plant species", "Animal species")),
          legend.title=" ") 
# figS3 %>% ggsave(filename=here("figures/corr_maxdispersalrange_rggcon.pdf"),width=8,height=8)

#### Figure S4: Food webs before (a+c) and after invasion (b+d) ###
plot_web <- function (A  = web, vertexSizeFactor = 3, title = webname) {
  set.seed(54321)
  A <- ifelse(A[,]>0, 1, 0) %>% unname
  g <- graph_from_adjacency_matrix(t(A))
  deg <- degree(g, mode="all") #mass #disp
  V(g)$frame.color <- "white"
  V(g)$color <- "orange"
  
  tl <- TrophInd(get.adjacency(g, sparse = F))
  V(g)$label <- round(tl$TL)
  
  lMat <- matrix(nrow = vcount(g), ncol = 2)
  lMat[, 2] <- tl$TL
  lMat[, 1] <- runif(vcount(g))
  
  colourCount <- tl$TL %>% round %>% unique %>% length
  getPalette  <- colorRampPalette(brewer.pal((11),"RdYlBu"))
  cpalette <- rev(getPalette(colourCount+1))
  cpalette <- c("springgreen4", cpalette[-4])
  V(g)$color <- cpalette[round(tl$TL)]
  V(g)$label <- round(tl$TL)
  plot(g, edge.width = 0.3, edge.arrow.size = 0.2, vertex.label.color = "white", 
       edge.color = "grey30", edge.curved = 0.2, layout = lMat, main = title)
  }

webfiles <- c(list.files(path=here("data/output/"), pattern = "^web_41516_1_1_1000_0.csv$", all.files = TRUE,
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE), 
              list.files(path=here("data/output/"), pattern = "^web_88516_1_20_1000_0.csv$", all.files = TRUE,
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE), 
              list.files(path=here("data/webs"), pattern = "^web_41516_99.csv$", all.files = TRUE,
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE), 
              list.files(path=here("data/webs"), pattern = "^web_88516_99.csv$", all.files = TRUE,
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE))

for(webfile in webfiles){
  pdfname <- webfile %>% basename %>% sub(".csv", ".pdf", .) %>% 
    paste0(here("figures/"), .)
  # pdf(pdfname)
  web <- webfile %>% read.table %>% as.matrix 
  if(pdfname %>% str_detect(.,pattern="web_41516")) webname <- "food web 1"
  if(pdfname %>% str_detect(.,pattern="web_88516")) webname <- "food web 2"
  if(pdfname %>% str_detect(.,pattern="web_7516")) webname <- "food web 3"
  if(pdfname %>% str_detect(.,pattern="web_8516")) webname <- "food web 4"
  if(pdfname %>% str_detect(.,pattern="web_42516")) webname <- "food web 5"
  mytitle <- ifelse(pdfname %>% str_detect(.,pattern="_99"), 
                    paste(webname, "after invasion"), 
                    paste(webname, "before invasion"))
  plot_web(web, vertexSizeFactor = 7, title = mytitle)
  # dev.off()
  }


#### Figure S5: Dispersal rate ~ net growth rate ####
#Function to calculate dispersal rate based on net growth rate
dispersal.rate <- function(G=0.9, a=0.1, b=10, c=0.141){
  nom <- a
  denom <- 1 + exp(b*(-c-G))
  D <- (nom/denom)
  return(D)
}

#Parameters
G <- seq(-1, 1, 0.01) #net growth rates
a <- 0.1 #maximum dispersal rate
b_b <- 10 #shape parameter plants (increasing S-fucntion)
b_c <- -10 #shape parameter consumers (decreasing S-function)

# pdf(here("figures/dispersal_rate.pdf"), width=16, height=8)

par(mfrow=c(1,2))
par(pty="s")

#Plants
mat = outer(X=G, Y=b_b, function(X,Y)  dispersal.rate(G=X,a=a,b=Y,c=0.138))
matplot(G,mat,type="l",ylab="Dispersal rate",
        xlab="Net growth rate",cex.lab=2.5,
        lwd = 2, col = spcols[1], lty=1, cex.axis = 1.5)
segments(-0.138,0,-0.138,0.05,xpd=T,lty=2,lwd=2,col=spcols[1])
segments(-1,0.05,-0.138,0.05,xpd=T,lty=3,lwd=2,col="darkgrey")

#Animals
mat = outer(X=G, Y=b_c, function(X,Y)  dispersal.rate(G=X,a=a,b=Y,c=0.141))
matplot(G,mat,type="l",ylab="",xlab="Net growth rate",cex.lab=2.5,
        lwd = 2, col = spcols[2], lty=1, cex.axis = 1.5)
segments(-0.141,0,-0.141,0.05,xpd=T,lty=2,lwd=2,col=spcols[2])
segments(-1,0.05,-0.141,0.05,xpd=T,lty=3,lwd=2,col="darkgrey")

# dev.off()

#### Figure S6: Emerged beta-diversity t = 10.000 ~ init. beta-diversity t = 0 ####
figS6 <- ggscatter(dat_plot, x = "beta", y = "beta.post",
          add = "none", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Initialized beta-diversity (t=0)", ylab = "Emerged beta-diversity (t=10.000)", 
          color="landscape.type", palette=c("blue","darkgreen"),
          # color=c("landscape.type","nutrient.supply"), 
          facet.by = c("landscape.type","nutrient.scenario"), 
          legend.title="spatial habitat configuration") 
# figS6 %>% ggsave(filename=here("figures/corr_beta.pdf"),width=16,height=8)

## Figure S7C: Species richness ~ landscape type + nutrient supply ###
x <- dat_plot %>% pull(landscape.type)
xu <- levels(x)
pal <- gg_color_hue(xu %>% length+1)[2:4]
y <- dat_plot %>% pull(gamma.post)

figS7c <- dat_plot %>%
  ggplot(aes(x=as.factor(x), y=y, alpha=0.4)) +
  stat_summary(aes(colour=as.factor(x), fill=as.factor(x)), alpha=0.4, 
               show.legend = FALSE, position = position_dodge(1), 
               fun.data = quantiles_95, geom="boxplot") + 
  scale_y_continuous(name = "Species richness\n") +
  geom_hline(aes(yintercept=21), alpha=0.3, colour="black", lwd=1, linetype="dashed") +
  scale_x_discrete(name = "") +
  scale_colour_manual(values=bg) +
  scale_fill_manual(values=bg) +
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
         col = guide_legend(reverse=FALSE)) + 
  stat_compare_means(comparisons = list(c("random", "clustered")), 
                     label = "p.signif", method = "t.test")
# figS7c %>% ggsave(filename = here("figures/figS7c.pdf"), width=16, height=12)

#### Figure S8: Species richness ~ landscape type + nutrient supply ###
datin <- dat_plot %>% drop_na(rem.spp)
datin$x <- datin %>% pull(rem.spp)
xu <- datin$x %>% unique %>% sort

figS8 <- datin %>% ggplot(aes(x=x, y=fraction.of.invaded.patches)) +  
  geom_boxplot(aes(group=x, colour=invspp.type, fill=invspp.type), alpha = 0.4, 
               show.legend = FALSE, position = position_dodge2(preserve = "single"), notch=FALSE,
               stat = "summary", fun.data = quantiles_95) +
  scale_y_continuous(name = "Fraction invaded patches\n",
                     breaks = seq(0, 100, 50),
                     limits=c(0, 100), 
                     labels = function(b) { paste0(round(b * 1, 0), "%")}) + 
  scale_x_discrete(name = "\nInvasive species", 
                   breaks = xu,
                   limits = factor(xu), 
                   labels = xu) +
  scale_colour_manual(values = spcols) +
  scale_fill_manual(values = spcols) +
  facet_grid(landscape.type ~ nutrient.scenario, 
             labeller=labeller(.rows = label_value, .cols = label_value)) + 
  ggpubr::theme_pubr() +
  theme(text = element_text(size=20), 
        strip.text.y = element_text(size=40),
        strip.text.x = element_text(size=40),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        axis.text.x = element_text(size = 14, angle=90),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 44),
        axis.title.x = element_text(size = 44))

# figS8 <- ggsave(filename = here("figures/invspp.pdf"), width=16, height=12)

#### Figure S9A: Emerged log10 biomass density t = 10.000 ~ init. log10 biomass density t = 0 ####
figS9a <- ggscatter(dat_plot, x = "log10tott0", y = "log10totpost",
                     add = "reg.line", conf.int = TRUE,
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "Initialized biomass density [log10] (t=0)",
                     ylab = "Emerged biomass density [log10] (t=10.000)",
                     color="webno", fill="webno", facet.by = "webno",
                     legend.title = " ") 
# figS9a %>% ggsave(filename=here("figures/corr_biomasstot.pdf"),width=16,height=8)

#### Figure S9B: Emerged log10 invader biomass density t = 10.000 ~ init. log10 invader biomass density t = 5.000 ####
figS9b <- ggscatter(dat_plot, x = "log10invsppint", y = "log10invspppost",
                     add = "reg.line", conf.int = TRUE,
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "Initialized invader biomass density [log10] (t=5.000)",
                     ylab = "Emerged invader biomass density [log10] (t=10.000)",
                     color = "webno", fill="webno", facet.by="webno", 
                     legend.title = "") 
# figS9b %>% ggsave(filename=here("figures/corr_biomassinvspp.pdf"),width=16,height=8)
