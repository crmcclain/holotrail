########################required packages########################
#################################################################

    patch_packages <- c("dplyr", "ggplot2", "vegan", "readxl","gridExtra", "cluster", 
                        "factoextra", "tibble", "ggrepel", "ggdendro", "betapart", 
                        "modelsummary", "factoextra", "scatterpie")
    lapply(patch_packages, require, character.only = TRUE)
#############################themes##############################
#################################################################
    theme_craig <- function () { 
      theme_bw(base_size=12, base_family="Helvetica") %+replace% 
        theme(
          # change stuff here
          axis.line = element_line(colour = "darkgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(face="bold"),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0),
          legend.position="none")
    }
    
    theme_craig_legend <- function () { 
      theme_bw(base_size=12, base_family="Helvetica") %+replace% 
        theme(
          # change stuff here
          axis.line = element_line(colour = "darkgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(face="bold"),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
    }

    
    study_color_pallette <- scale_color_manual(values = c("trail" = "#009E73",
                                                          "background" = "#0072B2"))
    
    study_color_pallette2 <- scale_fill_manual(values = c("trail" = "#009E73",
                                                          "background" = "#0072B2"))
    
    
    
############################load data############################
#################################################################  
    setwd("~/Dropbox/Return of the Woodfall/Projects/Poop Trail")
    
    
    corebytaxa<- read.csv("poop_trail.csv") %>%
      replace(is.na(.), 0)
    
    NCOL <- ncol(corebytaxa)
    spe <- corebytaxa[,3:NCOL] 
    
    coresummary<-corebytaxa %>%
      select(Core, trail)

############################PCOA############################
#################################################################  

    
    
    # Ordination  (PCoA)
    spe.norm <- decostand(spe, "normalize")
    spe.euclid <- vegdist (spe.norm, "euc")
    spe.pcoa <- cmdscale(spe.euclid, eig=TRUE)
    
    
    
    #get dim
    coresummary$DIM1 <-scores(spe.pcoa)[,1]
    coresummary$DIM2 <-scores(spe.pcoa)[,2]   
    
    spe.pcoa$eig[1]/sum(spe.pcoa$eig)
    spe.pcoa$eig[2]/sum(spe.pcoa$eig)
    
    
    #get species scores
    species.scores <- data.frame(wascores(spe.pcoa$points[,1:2], spe))
    #limit to upper quantiles
    species.scores2 <- species.scores %>%
      filter(
        X1<= quantile (species.scores$X1, probs=.05, na.rm=TRUE) | 
          X1>= quantile (species.scores$X1, probs=.95,na.rm=TRUE) | 
          X2<= quantile (species.scores$X2, probs=.05,na.rm=TRUE) | 
          X2>= quantile (species.scores$X2, probs=.95,na.rm=TRUE) )  %>%
      rownames_to_column('name')
    
    #And Finally the Plot Together
    pc1 <- ggplot(data=coresummary, aes(DIM1, DIM2))+
      geom_point(alpha=0.9, pch=21, cex=5, aes(fill=trail))+
      #geom_text(aes(label=Core, color=trail))+
      coord_equal()+
      theme_craig()+
      #geom_text_repel(data=species.scores2, aes(x=X1, y=X2, label=name), color="red")+
      study_color_pallette+
      study_color_pallette2
    
    pc1


    ###adonis
    spe.h <- decostand (spe, "hellinger")
    adonis2(spe.h ~ trail, data = coresummary)  

############################Clustering############################
#################################################################  
    
    #Ward’s Minimum Variance Clustering
    spe.ch.ward <- hclust(spe.euclid, method="ward.D")
    spe.ch.ward$height <- sqrt(spe.ch.ward$height)
    spe.ch.ward$labels <- coresummary$trail
    spe.ch.ward$labels <- coresummary$Core

    plot(spe.ch.ward)
    
    
############################LCBD############################
#################################################################  
    #within group beta
    require(adespatial)
    

    
    overall_beta <-beta.div(spe, method="hellinger" ,nperm = 999)
    overall_beta
    coresummary$LCBD<- overall_beta$LCBD
    
    ggplot(data=coresummary, aes(y=LCBD, x=trail, group=trail, fill=trail))+
      geom_boxplot()+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2
    
    summary(aov(data=coresummary, LCBD~trail))
    
#######################Calculate Diversity########################
#################################################################     
    
    
    
    #diversity time
    
    corebytaxadiv <-   corebytaxa %>% 
      mutate(Abundance = rowSums( corebytaxa[,3:NCOL]), #Abundance                
             S = specnumber(corebytaxa[,3:NCOL]), #Richness 
             H = diversity( corebytaxa[,3:NCOL],index="shannon"), #Shannon's H
             Simp = diversity( corebytaxa[,3:NCOL],index="simpson"), #Simpson J
             log10Abundance=log10(Abundance)) %>% #log Abundance
      select(Abundance, log10Abundance, S, H, Simp,Core)
    
    #joining diversity sheets together
    
    coresummary <- left_join(coresummary, corebytaxadiv, by=c("Core"="Core")) 
 
    
    
    div1 <-   ggplot(data=coresummary, aes(y=S, x=trail, group=trail, fill=trail))+
      geom_boxplot()+
      study_color_pallette+
      study_color_pallette2+
      theme_craig()+
      ylab("Species Richness (S)")
    
    div2 <-   ggplot(data=coresummary, aes(y=H,  x=trail, group=trail, fill=trail))+
      geom_boxplot()+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Shannon-Weiner H'")
    
    div3 <-   ggplot(data=coresummary, aes(y=Simp,  x=trail, group=trail, fill=trail))+
      geom_boxplot()+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Simpson J'")
    
    div4 <-   ggplot(data=coresummary, aes(y=log10Abundance,  x=trail, group=trail, fill=trail))+
      geom_boxplot()+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Log 10 Abundnace")
    
    grid.arrange(div1, div2, div3, div4, nrow=2)
    
    summary(aov(data=coresummary, S~trail))
    t.test(data=coresummary, S~trail,  alternative = c("greater"))
    summary(aov(data=coresummary, H~trail))
    t.test(data=coresummary, H~trail,  alternative = c("greater"))
    summary(aov(data=coresummary, Simp~trail))
    t.test(data=coresummary, Simp~trail,  alternative = c("greater"))
    summary(aov(data=coresummary, log10Abundance~trail))
    t.test(data=coresummary, log10Abundance~trail,  alternative = c("greater"))
    
#######################indval########################
#################################################################   
    require(labdsv)
    (iva <- indval(spe, coresummary$trail)) 
    
    # Table of the significant indicator species
    gr <- iva$maxcls[iva$pval <= 0.1]
    iv <- iva$indcls[iva$pval <= 0.1]
    pv <- iva$pval[iva$pval <= 0.1]
    
    rf1 <- iva$relfrq[iva$pval <= 0.1,1]
    rf2 <- iva$relfrq[iva$pval <= 0.1,2]
    
    ra1 <- iva$relabu[iva$pval <= 0.1,1]
    ra2 <- iva$relabu[iva$pval <= 0.1,2]
    
    fr <- apply(spe > 0, 2, sum)[iva$pval <= 0.1]
    fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr, relfrq1=rf1, relfrq2=rf2, relabu1=ra1, relabu2=ra2)
    fidg <- fidg[order(fidg$group, -fidg$indval),]
    fidg
    
    