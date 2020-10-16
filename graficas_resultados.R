multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

if(!require('ggplot2')) install.packages('ggplot2')
require('ggplot2')

#install.packages("hrbrthemes")

# Libraries
library(tidyverse)
library(hrbrthemes)
#library(kableExtra)
#options(knitr.table.format = "html")
#library(babynames)
#library(streamgraph)
library(viridis)
library(DT)
#library(plotly)

data<-read.table("C:/Users/Carmen C/Documents/R/Proyecto-Maestria/Tablas/resultados.txt", 
                header = TRUE)

head(data)


# Plot de las cartas
# Filter data
data_g <- data %>%
  filter(carta %in% c("VP","VSSIWL", "VSSI", "VSI", "VSS", "X")) %>%
  filter(phi==0.2)%>%
  filter(psi==0.5)%>%
  filter(delta!=0.1)
p1<-(ggplot(data=data_g, aes(x=delta, y=ATS1, group=carta, color=carta)) +
       geom_point() + geom_line() +
       scale_color_manual(values = c("red", "blue", "green", "cyan", "darkgoldenrod", "darkorchid")) +
       theme(legend.position="none",
             plot.title = element_text(size=14)) +
       ggtitle("ATS1 para phi=0.2") +
       theme_ipsum())

data_g <- data %>%
  filter(carta %in% c("VP","VSSIWL", "VSSI", "VSI", "VSS", "X")) %>%
  filter(phi==0.4)%>%
  filter(psi==0.5)%>%
  filter(delta!=0.1)
p2<-(ggplot(data=data_g, aes(x=delta, y=ATS1, group=carta, color=carta)) +
       geom_point() + geom_line() +
       scale_color_manual(values = c("red", "blue", "green", "cyan", "darkgoldenrod", "darkorchid")) +
       theme(legend.position="none",
             plot.title = element_text(size=14)) +
       ggtitle("ATS1 para phi=0.4") +
       theme_ipsum())

data_g <- data %>%
  filter(carta %in% c("VP","VSSIWL", "VSSI", "VSI", "VSS", "X")) %>%
  filter(phi==0.6)%>%
  filter(psi==0.5)%>%
  filter(delta!=0.1)
p3<-(ggplot(data=data_g, aes(x=delta, y=ATS1, group=carta, color=carta)) +
       geom_point() + geom_line() +
       scale_color_manual(values = c("red", "blue", "green", "cyan", "darkgoldenrod", "darkorchid")) +
       theme(legend.position="none",
             plot.title = element_text(size=14)) +
       ggtitle("ATS1 para phi=0.6") +
       theme_ipsum())

data_g <- data %>%
  filter(carta %in% c("VP","VSSIWL", "VSSI", "VSI", "VSS", "X")) %>%
  filter(phi==0.8)%>%
  filter(psi==0.5)%>%
  filter(delta!=0.1)
p4<-(ggplot(data=data_g, aes(x=delta, y=ATS1, group=carta, color=carta)) +
       geom_point() + geom_line() +
       scale_color_manual(values = c("red", "blue", "green", "cyan", "darkgoldenrod", "darkorchid")) +
       theme(legend.position="none",
             plot.title = element_text(size=14)) +
       ggtitle("ATS1 para phi=0.8") +
       theme_ipsum())

multiplot(p1, p2, p3, p4, cols = 2)

# Plot VP
data_VP <- data %>%
  filter(carta %in% ("VP")) %>%
    filter(psi==0.5)%>%
      #filter(phi!=0.99)%>%
        filter(delta!=0.1)

data_VP

ggplot(data=data_VP, aes(x=delta, y=ATS1, group=phi, color=as.factor(phi))) +
  geom_point() + geom_line() +
scale_color_manual(values = c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", 
                              "orange")) +
  theme(legend.position="none",
        plot.title = element_text(size=14)) +
  ggtitle("ATS1 de la carta sintética VP para psi = 0.5") +
  theme_ipsum()
