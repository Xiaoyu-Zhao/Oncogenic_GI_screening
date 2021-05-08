library(cowplot)
library(ggplot2)
gseaPlot = function(gsdata, geneSetID, title = "", color, base_size = 26, dotsize,
                    rel_heights=c(1.5, .5, 1), subplots = 1:3, pvalue_table = FALSE, ES_geom="line",
                    y_min, y_max, legend.position) {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
 p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size)
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description), size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description), size=dotsize, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer 
  p.res <- p.res + ylab("RES" )+
    theme(panel.border = element_blank(), 
          panel.grid = element_blank(),
          panel.background = element_blank(),
          #panel.grid= element_line(linetype = "dashed"),
          axis.line = element_line(colour = "black"),
          legend.position = 'none',
          legend.title = element_blank(),
          #legend.title = element_text(size = 20, face= "bold"),
          #legend.text = element_text(size=24, color = "black"),
          axis.text.y  = element_text(size=38, color = "black"),
          axis.text.x =element_text(size=38, color = "black"),
          axis.title=element_text(size=40, color = "black"),
          plot.title = element_text(size=40,color = "black"),
          plot.margin=margin(t = .5, r = .5, b=.5, l=.5, unit="cm")
          #strip.text.y  = element_text(size = 26, face="bold")
          ##strip.text is for facet
    ) +scale_x_continuous( limits = c(0,16000), breaks = seq(0, 15000, by = 5000)) +
    scale_y_continuous(limits = c(y_min,y_max), breaks = round(seq(y_min, y_max, by = 0.2),1))+
    expand_limits(x = 0, y = 0)
  
  i <- 0
  terms = unique(gsdata$Description)
  terms = sort(terms[!is.na(terms)], decreasing = TRUE, na.last = TRUE)
  for (term in terms) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1.2
  }
  p2 <- ggplot(gsdata, aes_(x = ~x))+
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax), color = "black")+ 
    xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) +
    theme(panel.border = element_blank(), 
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)+ 
    theme(plot.title = element_text(size=36,hjust=0.5, vjust = 2))
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
   
  }
  
  plotlist <- list(p.res, p2)
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(plot.margin=margin(t=.5, r = .5, b=.5, l=.5, unit="cm"))
  
  if (length(subplots) == 1)
    {return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2, l=.2, unit="cm")))}
  
  if (length(rel_heights) > length(subplots))
    {rel_heights <- rel_heights[subplots]}
  
  plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights = rel_heights)
}
