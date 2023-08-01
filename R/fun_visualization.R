
## Visualization helpers

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
getPaletteRdBu <- colorRampPalette(brewer.pal(8, "RdBu"))
getPalette_uchi <- colorRampPalette(c(ggsci::pal_uchicago(palette = "light")(9)))
getPalette_uchi2 <- colorRampPalette(c(ggsci::pal_uchicago(palette = "default")(9)))
getPalette_jco   <- colorRampPalette(c(ggsci::pal_jco()(10)))
getPaletteCustom <- function(n) {
  palette <- c("dodgerblue1", "skyblue4", "chocolate1", "seagreen4",
               "bisque3", "red4", "purple4", "mediumpurple3",
               "maroon", "dodgerblue4", "skyblue2", "darkcyan",
               "darkslategray3", "lightgreen", "bisque",
               "palevioletred1", "black", "gray79", "lightsalmon4",
               "darkgoldenrod1")
  if (n > length(palette))
    warning('palette has duplicated colours')
  rep(palette, length.out=n)
}

facets_nice2 <- theme(strip.text.x = element_text(size=15, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))
facets_nice3 <- theme(strip.text.x = element_text(size=8, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))
facets_nice4 <- theme(strip.text.x = element_text(size=11, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))

add_block <- theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line=element_blank())

plot_feature_nice <- function(){labs(x = "", y = "") & 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) & 
  add_block & theme(legend.position = "none")}

plot_feature_nice <- theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
