
#Figure 1b

setwd("~/Desktop/Project_ITS_CNV")
library(data.table)
library(plotrix)

#patch to fix gap.boxplot
gap.boxplot<-function (x,...,gap=list(top=c(NA,NA),bottom=c(NA,NA)),
                       range = 1.5, width = NULL, varwidth = FALSE, notch = FALSE, 
                       outline = TRUE, names, xlim = NA, ylim = NA, plot = TRUE, 
                       border = par("fg"), col = NULL, log = "", axis.labels = NULL, 
                       axes = TRUE, pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5), 
                       horizontal = FALSE, add = FALSE, at = NULL, main = NULL) {
  
  if (!is.na(gap$top[1])) 
    if (gap$top[1] > gap$top[2]) 
      gap$top <- rev(gap$top)
    if (!is.na(gap$bottom[1])) 
      if (gap$bottom[1] > gap$bottom[2]) 
        gap$bottom <- rev(gap$bottom)
      if (is.na(ylim[1])) {
        bxpt <- boxplot(x, ..., range = range, plot = FALSE)
        ylim <- range(c(bxpt$stats, bxpt$out))
      }
      else bxpt <- boxplot(x, ..., ylim = ylim, range = range,plot = FALSE)
      bxgap <- bxpt
      if (!is.na(gap$top[1])) {
        bxgap$stats[bxgap$stats > gap$top[1] & bxgap$stats < 
                      gap$top[2]] <- NA
        if (any(is.na(bxgap$stats))) 
          stop("gap cannot include the median, interquartiles or the staples")
        topdiff <- diff(gap$top)
        bxgap$stats[bxgap$stats > gap$top[2]] <- bxgap$stats[bxgap$stats > 
                                                               gap$top[2]] - topdiff
        intopgap <- bxgap$out > gap$top[1] & bxgap$out < gap$top[2]
        bxgap$out[intopgap] <- NA
        abovetop <- which(bxgap$out > gap$top[2])
        bxgap$out[abovetop] <- bxgap$out[abovetop] - topdiff
        rangetop <- gap$top[1]
        ylim[2] <- ylim[2] - topdiff
      }
      else rangetop <- ylim[2]
      if (!is.na(gap$bottom[1])) {
        bxgap$stats[bxgap$stats > gap$bottom[1] & bxgap$stats < 
                      gap$bottom[2]] <- NA
        if (any(is.na(bxgap$stats))) 
          stop("gap cannot include the median, interquartiles or the staples")
        bottomdiff <- diff(gap$bottom)
        bxgap$stats[bxgap$stats < gap$bottom[1]] <- bxgap$stats[bxgap$stats < 
                                                                  gap$bottom[1]] + bottomdiff
        bxgap$out[bxgap$out > gap$bottom[1] & bxgap$out < gap$bottom[2]] <- NA
        belowbottom <- which(bxgap$out < gap$bottom[1])
        bxgap$out[belowbottom] <- bxgap$out[belowbottom] + bottomdiff
        rangebottom <- gap$bottom[2]
        ylim[1] <- ylim[1] + bottomdiff
      }
      else rangebottom <- ylim[1]
      if (any(is.na(bxgap$out))) warning("At least one outlier falls into a gap")
      nboxes <- dim(bxgap$stats)[2]
      if (is.na(xlim[1])) {
        xlim <- c(0.5, nboxes + 0.5)
        at <- 1:nboxes
      }
      bxgap$group <- at
      plot(0,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab="",ylab="",main=main)
      plotlim <- par("usr")
      box()
      if (axes) axis(1, labels = bxpt$names, at = at)
      midticks <- pretty(c(rangebottom, rangetop))
      if(axes) axis(2,at=midticks[midticks > rangebottom & midticks < rangetop])
      if (is.null(width)) width <- pars$boxwex
      rect(at - width/2, bxgap$stats[2, ], at + width/2, bxgap$stats[4, 
                                                                     ], border = border, col = col)
      if (notch) {
        ymult <- getYmult()
        if (is.null(col)) 
          boxcol <- "white"
        else boxcol <- col
        rect(at - width/1.95, bxgap$conf[1, ], at + width/1.95, 
             bxgap$conf[2, ], border = NA, col = boxcol)
        insets <- (bxgap$conf[2, ] - bxgap$conf[1, ]) * pars$boxwex/ymult
        median.left <- ((at - width/2) + insets)
        median.right <- ((at + width/2) - insets)
        segments(at - width/2, bxgap$conf[1, ], median.left, 
                 bxgap$stats[3, ], col = border)
        segments(at - width/2, bxgap$conf[2, ], median.left, 
                 bxgap$stats[3, ], col = border)
        segments(median.right, bxgap$stats[3, ], at + width/2, 
                 bxgap$conf[1, ], col = border)
        segments(median.right, bxgap$stats[3, ], at + width/2, 
                 bxgap$conf[2, ], col = border)
      }
      else {
        median.left <- at - width/2
        median.right <- at + width/2
      }
      segments(median.left, bxgap$stats[3, ], median.right, bxgap$stats[3, 
                                                                        ], lwd = 2, col = border)
      segments(at, bxgap$stats[1, ], at, bxgap$stats[2, ], lty = 2, 
               col = border)
      segments(at, bxgap$stats[4, ], at, bxgap$stats[5, ], lty = 2, 
               col = border)
      segments(at - pars$staplewex * width/2, bxgap$stats[1, ], 
               at + pars$staplewex * width/2, bxgap$stats[1, ], col = border)
      segments(at - pars$staplewex * width/2, bxgap$stats[5, ], 
               at + pars$staplewex * width/2, bxgap$stats[5, ], col = border)
      if (!is.na(gap$top[1])) topadjust <- diff(gap$top)
      else topadjust <- 0
      if (!is.na(gap$bottom[1])) bottomadjust <- diff(gap$bottom)
      else bottomadjust <- 0
      if (!is.null(axis.labels)) 
        axis(2, labels = axis.labels, at = c(axis.labels[1] + 
                                               bottomadjust, axis.labels[2] - topadjust))
      if (!is.na(gap$top[1])) axis.break(2, gap$top[1], style = "gap")
      if (!is.na(gap$bottom[1])) 
        axis.break(2, gap$bottom[2] - diff(plotlim[3:4]) * 0.02, 
                   style = "gap")
      print(bxgap)
      points(bxpt$group,bxgap$out)
      invisible(bxgap)
}


rDNA_by_taxa<- read.csv("rDNA_by_taxa.csv", sep = ",", header = TRUE)
#rDNA_by_taxa <- rDNA_by_taxa_full[ rDNA_by_taxa_full$Genus!="Basidiobolus",]


#smaller df of only species with multiple representatives 
#get range for each group
SE_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$SE) | duplicated(rDNA_by_taxa$SE, fromLast = TRUE),]
min_by_SE<-aggregate(CN ~ SE, SE_table, function(x) min(x))
max_by_SE<-aggregate(CN ~ SE, SE_table, function(x) max(x))
SE_range<- cbind(min=min_by_SE$CN, max=max_by_SE$CN)
dif<- abs(SE_range[,1] - SE_range[,2])
SE_range<- cbind(SE_range, dif)
rownames(SE_range) <- min_by_SE$SE
#get median of range values 
SE_difference<- median(abs(SE_range[,1] - SE_range[,2]))

#smaller df of only genera with multiple representatives 
#get range for each group 
Genus_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$Genus) | duplicated(rDNA_by_taxa$Genus, fromLast = TRUE),]
min_by_Genus<-aggregate(CN ~ Genus, Genus_table, function(x) min(x))
max_by_Genus<-aggregate(CN ~ Genus, Genus_table, function(x) max(x))
Genus_range<- cbind(min=min_by_Genus$CN, max=max_by_Genus$CN)
dif<- abs(Genus_range[,1] - Genus_range[,2])
Genus_range<- cbind(Genus_range, dif)
rownames(Genus_range) <- min_by_Genus$Genus
#get median of range values 
Genus_difference<- median(abs(Genus_range[,1] - Genus_range[,2]))


#smaller df of only families with multiple representatives 
#get range for each group
Family_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$Family) | duplicated(rDNA_by_taxa$Family, fromLast = TRUE),]
min_by_Family<-aggregate(CN ~ Family, Family_table, function(x) min(x))
max_by_Family<-aggregate(CN ~ Family, Family_table, function(x) max(x))
Family_range<- cbind(min=min_by_Family$CN, max=max_by_Family$CN)
dif<- abs(Family_range[,1] - Family_range[,2])
Family_range<- cbind(Family_range, dif)
rownames(Family_range) <- min_by_Family$Family
#get median of range values 
Family_difference<- median(abs(Family_range[,1] - Family_range[,2]))

#smaller df of only orders with multiple representatives 
#get range for each group
Order_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$Order) | duplicated(rDNA_by_taxa$Order, fromLast = TRUE),]
min_by_Order<-aggregate(CN ~ Order, Order_table, function(x) min(x))
max_by_Order<-aggregate(CN ~ Order, Order_table, function(x) max(x))
Order_range<- cbind(min=min_by_Order$CN, max=max_by_Order$CN)
dif<- abs(Order_range[,1] - Order_range[,2])
Order_range<- cbind(Order_range, dif)
rownames(Order_range) <- min_by_Order$Order
#get median of range values 
Order_difference<- median(abs(Order_range[,1] - Order_range[,2]))

#smaller df of only classes with multiple representatives 
#get range for each group
Class_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$Class) | duplicated(rDNA_by_taxa$Class, fromLast = TRUE),]
min_by_Class<-aggregate(CN ~ Class, Class_table, function(x) min(x))
max_by_Class<-aggregate(CN ~ Class, Class_table, function(x) max(x))
Class_range<- cbind(min=min_by_Class$CN, max=max_by_Class$CN)
dif<- abs(Class_range[,1] - Class_range[,2])
Class_range<- cbind(Class_range, dif)
rownames(Class_range) <- min_by_Class$Class
#get median of range values 
Class_difference<- median(abs(Class_range[,1] - Class_range[,2]))

#smaller df of only phyla with multiple representatives 
#get range for each group
Phylum_table<- rDNA_by_taxa[duplicated(rDNA_by_taxa$Phylum) | duplicated(rDNA_by_taxa$Phylum, fromLast = TRUE),]
min_by_Phylum<-aggregate(CN ~ Phylum, Phylum_table, function(x) min(x))
max_by_Phylum<-aggregate(CN ~ Phylum, Phylum_table, function(x) max(x))
Phylum_range<- cbind(min=min_by_Phylum$CN, max=max_by_Phylum$CN)
dif<- abs(Phylum_range[,1] - Phylum_range[,2])
Phylum_range<- cbind(Phylum_range, dif)
rownames(Phylum_range) <- min_by_Phylum$Phylum
#get median of range values 
Phylum_difference<- median(abs(Phylum_range[,1] - Phylum_range[,2]))

#By mean
SE_difference.mean<- mean(abs(SE_range[,1] - SE_range[,2]))
Genus_difference.mean<- mean(abs(Genus_range[,1] - Genus_range[,2]))
Family_difference.mean<- mean(abs(Family_range[,1] - Family_range[,2]))
Order_difference.mean<- mean(abs(Order_range[,1] - Order_range[,2]))
Class_difference.mean<- mean(abs(Class_range[,1] - Class_range[,2]))
Phylum_difference.mean<- mean(abs(Phylum_range[,1] - Phylum_range[,2]))

#figure 1.B, without a brake 
boxplot(SE_range[,3], 
        Genus_range[,3], 
        Family_range[,3], 
        Order_range[,3], 
        Class_range[,3], 
        Phylum_range[,3], 
        range = 3, 
        outline=TRUE, 
        names = c("SE", "G", "F", "O", "C", "P"))
axis.break(2, from, breakcol="white", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")


# but we need to brake the axis, because that outlier makes it dificult to see the real data.
#this moves all the outliers into the wrong group??
from <- 300
to <- 1200
gap.boxplot(SE_range[,3], 
        Genus_range[,3], 
        Family_range[,3], 
        Order_range[,3], 
        Class_range[,3], 
        Phylum_range[,3], 
        gap=list(top=c(from,to),bottom=c(NA,NA)),
        range = 3, 
        outline = TRUE,
        names = c("SE", "G", "F", "O", "C", "P"))
axis.break(2, from, breakcol="white", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2,at=c(1,50,150,200,350,400, 450),labels=c("1","50","150","200","1200","1250", "1300"))

#if we change range = 3 to range 0, the wiskers extend to the correct group outliers
#.. why does this brake when we need both wiskers and outliers to be explicit.
from <- 300
to <- 1200
plotrix::gap.boxplot(SE_range[,3], 
            Genus_range[,3], 
            Family_range[,3], 
            Order_range[,3], 
            Class_range[,3], 
            Phylum_range[,3], 
            gap=list(top=c(from,to),bottom=c(NA,NA)),
            range = 3, 
            outline = TRUE,
            names = c("SE", "G", "F", "O", "C", "P"))
axis.break(2, from, breakcol="white", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2,at=c(1,50,150,200,350,400, 450),labels=c("1","50","150","200","1200","1250", "1300"))


