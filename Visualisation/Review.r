Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/Cowplot.reviex/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
x <- get(load(file = paste0(Allmod[i],".rda")))
dput(list.files(pattern=".rds"))

NAME <- c("ABIALB","BETPEN","CASSAT","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEROB")

for (i in c(1:10)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}

library(cowplot)
library(ggplot2)

ptest <- pABIALB+theme(legend.position = c(0.5,0.55),legend.key.size = unit(2.5,"cm"),legend.text = element_text(size=18),legend.title = element_text(size=20))+
  scale_color_manual("Drought index effect according to Marginality levels:",values = c("black", "red", "blue"),labels=(c("Core","Trailing Edge","Leading Edge")))+ #change legend when one cat is removed
  scale_shape_manual("Drought index effect according to Marginality levels:",values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))+
  guides(color = guide_legend(title.position = "top",title.hjust=0.5))

#scale_shape_manual("Couille", values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))
legend <- get_legend(ptest)
plot(legend)

# part 1 occurrence 

pall <- plot_grid(pABIALB+theme(legend.position = "none",
              # axis.text.x=element_blank(),
              # axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              # axis.title.y=element_blank(),
              plot.title = element_blank(),
              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              
              pPICABI+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              pPINSYL+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              pCASSAT+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                             axis.title.x=element_blank(),
                            # axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              pPINPIN+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              pPINNIG+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                             axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              pPOPTRE+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                             axis.title.x=element_blank(),
                            # axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.8,0.2,0.5,0),"cm")),
              pQUEROB+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            # axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.8,0.2,0.5,0),"cm")),
              pBETPEN+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            # axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.8,0.2,0.5,0),"cm")),
              pPINHAL+theme(legend.position = "none",
                            # axis.text.x=element_blank(),
                            # axis.text.y=element_blank(),
                            # axis.title.x=element_blank(),
                            # axis.title.y=element_blank(),
                            plot.title = element_blank(),
                            plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
              legend,NULL,nrow = 4,ncol = 3,label_size = 16,labels = c("a) ABIALB",
                                                       "b) PICABI",
                                                       "c) PINSYL",
                                                       "d) CASSAT",
                                                       "e) PINPIN","f) PINNIG","g) POPTRE","h) QUEROB", "i) BETPEN", "j) PINHAL"),label_x = 0,label_y = 1.02)
pall   
pall2 <- pall + draw_line(x = c(0,1),y = c(0.505,0.505),color = "black", size = 1,linetype="dashed")
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Fig3_S5.ReviewV2.png"),plot = pall2,base_width = 6, base_height = 4, dpi = 300,units = "in",nrow=4,ncol=3)



p1 <- plot_grid(pABIALB+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                axis.title.x=element_blank(),
                                # axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  
                  pPICABI+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  pPINSYL+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  pCASSAT+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                #axis.title.x=element_blank(),
                                # axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  pPINPIN+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                #axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  pPINNIG+coord_cartesian(ylim = c(0,0.85))+theme(legend.position = "none",
                                # axis.text.x=element_blank(),
                                # axis.text.y=element_blank(),
                                #axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                plot.title = element_blank(),
                                plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                  nrow = 2,ncol = 3,label_size = 16,
                  labels = c("a) ABIALB","b) PICABI","c) PINSYL","d) CASSAT","e) PINPIN","f) PINNIG"),label_x = 0,label_y = 1.025)
                             
                             
  
p2 <- plot_grid(NULL,pPOPTRE+coord_cartesian(ylim = c(0,0.55))+theme(legend.position = "none",
                                           # axis.text.x=element_blank(),
                                           # axis.text.y=element_blank(),
                                           axis.title.x=element_blank(),
                                           # axis.title.y=element_blank(),
                                           plot.title = element_blank(),
                                           plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                             pQUEROB+coord_cartesian(ylim = c(0,0.55))+theme(legend.position = "none",
                                           #axis.text.x=element_blank(),
                                           # axis.text.y=element_blank(),
                                           axis.title.x=element_blank(),
                                           axis.title.y=element_blank(),
                                           plot.title = element_blank(),
                                           plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),NULL,NULL,
                             pBETPEN+coord_cartesian(ylim = c(0,0.55))+theme(legend.position = "none",
                                           # axis.text.x=element_blank(),
                                           # axis.text.y=element_blank(),
                                           # axis.title.x=element_blank(),
                                           #axis.title.y=element_blank(),
                                           plot.title = element_blank(),
                                           plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                             pPINHAL+coord_cartesian(ylim = c(0,0.55))+theme(legend.position = "none",
                                           # axis.text.x=element_blank(),
                                            #axis.text.y=element_blank(),
                                           # axis.title.x=element_blank(),
                                            axis.title.y=element_blank(),
                                           plot.title = element_blank(),
                                           plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),NULL,
                             nrow = 2,ncol = 4,label_size = 16,
                             labels = c("","g) POPTRE","h) QUEROB","","","i) BETPEN", "j) PINHAL",""),
                             label_x = 0,label_y = 1.025,rel_widths=c(0.5,1,1,0.5,0.5,1,1,0.5))
p2

p3 <- plot_grid(p1,legend,p2,nrow=3,ncol=1,rel_heights = c(0.99,0.25,0.9),scale = c(0.99,0.5,1))
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Fig3_S5.ReviewV5.png"),plot = p3,base_width = 6, base_height = 4.45, dpi = 1000,units = "in",nrow=4,ncol=3)


# first align the top-row plot (p3) with the left-most plot of the
# bottom row (p1)
plots <- align_plots(p1,p2, align = 'h', axis = 'l',nrow=2)
# then build the bottom row
bottom_row <- plot_grid(plots[[2]])
bottom_row
# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
p3 <- plot_grid(p1,legend,plots[[2]],nrow=3,ncol=1,rel_heights = c(1,0.2,1),scale = c(1,0.5,1))






## Figure 2
Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/Fig2.Review/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
x <- get(load(file = paste0(Allmod[i],".rda")))
dput(list.files(pattern=".rds"))

NAME <- c("ABIALB","FRAEXC","PINHAL","QUEILE","QUEPET","QUEROB")

for (i in c(1:10)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}

library(cowplot)
library(ggplot2)

ptest <- pABIALB+theme(legend.position = c(0.5,0.55),legend.key.size = unit(2.5,"cm"),legend.text = element_text(size=18),legend.title = element_text(size=20))+
  scale_color_manual("Temperature-related variable effect according to Marginality levels:",values = c("black", "red", "blue"),labels=(c("Core","Trailing Edge","Leading Edge")))+ #change legend when one cat is removed
  scale_shape_manual("Temperature-related variable effect according to Marginality levels:",values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))+
  guides(color = guide_legend(title.position = "top",title.hjust=0.5))

#scale_shape_manual("Couille", values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))
legend <- get_legend(ptest)
plot(legend)


p1 <- plot_grid(pABIALB+coord_cartesian(ylim = c(0,125))+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              # axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                
                pFRAEXC+coord_cartesian(ylim = c(0,125))+theme(legend.position = "none",
                              #axis.text.x=element_blank(),
                              #axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                pPINHAL+coord_cartesian(ylim = c(0,125))+xlim(25,31)+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              #axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                pQUEILE+coord_cartesian(ylim = c(0,125))+xlim(7.5,15)+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                               axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                pQUEPET+coord_cartesian(ylim = c(0,125))+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              #axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                pQUEROB+coord_cartesian(ylim = c(0,125))+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.2,0.5,0),"cm")),
                nrow = 3,ncol = 2,label_size = 16,
                labels = c("a) ABIALB","b) FRAEXC","c) PINHAL","d) QUEILE","e) QUEPET","f) QUEROB"),label_x = 0,label_y = 1.025)
p1

p3 <- plot_grid(p1,legend,nrow=2,ncol=1,rel_heights = c(0.99,0.15),scale = c(0.99,0.3))
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Fig4_S11.ReviewV1.pdf"),plot = p3,base_width = 4, base_height = 4, dpi = 600,units = "in",nrow=4,ncol=3)




## Fig S10

Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/Fig3.Review/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
x <- get(load(file = paste0(Allmod[i],".rda")))
dput(list.files(pattern=".rds"))

NAME <- c("BETPEN","FAGSYL","PICABI")

for (i in c(1:10)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}

library(cowplot)
library(ggplot2)

ptest <- pABIALB+theme(legend.position = c(0.5,0.55),legend.key.size = unit(2,"cm"),legend.text = element_text(size=16),legend.title = element_text(size=18))+
  scale_color_manual("Variable effect according to Marginality levels:",values = c("black", "red", "blue"),labels=(c("Core","Trailing Edge","Leading Edge")))+ #change legend when one cat is removed
  scale_shape_manual("Variable effect according to Marginality levels:",values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))+
  guides(color = guide_legend(title.position = "top",title.hjust=0.5))

#scale_shape_manual("Couille", values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))
legend <- get_legend(ptest)
plot(legend)


p1 <- plot_grid(pBETPEN+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              # axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                
                pFAGSYL+theme(legend.position = "none",
                              #axis.text.x=element_blank(),
                              #axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                pPICABI+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                nrow = 1,ncol = 3,label_size = 16,
                labels = c("a) BETPEN","b) FAGSYL","c) PICABI"),label_x = 0,label_y = 1.025)
p1
p3 <- plot_grid(p1,legend,nrow=2,ncol=1,rel_heights = c(0.95,0.4),scale = c(0.95,0.05))

cairo_pdf(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/TestFigS12.ReviewV1.pdf"),width = 16.5,height = 1.75*4,family="helvetica")
p3
dev.off()
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/FigS12.ReviewV1.pdf"),plot = p3,base_width = 5.5, base_height = 1.75, dpi = 300,units = "in",nrow=4,ncol=3)




### Figure suivante: 

## Fig S10

Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/S14/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
NAME <- c("FAGSYL","PINPINA","QUESUB")

for (i in c(1:3)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}

library(cowplot)
library(ggplot2)

ptest <- pFAGSYL+theme(legend.position = c(0.5,0.55),legend.key.size = unit(2,"cm"),legend.text = element_text(size=16),legend.title = element_text(size=18))+
  scale_color_manual("Variable effect according to Marginality levels:",values = c("black", "red", "blue"),labels=(c("Core","Trailing Edge","Leading Edge")))+ #change legend when one cat is removed
  scale_shape_manual("Variable effect according to Marginality levels:",values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))+
  guides(color = guide_legend(title.position = "top",title.hjust=0.5))

#scale_shape_manual("Couille", values=c(1,2,3),labels=(c("Core","Trailing Edge","Leading Edge")))
legend <- get_legend(ptest)
plot(legend)


p1 <- plot_grid(pFAGSYL+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              # axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                
                pPINPINA+theme(legend.position = "none",
                              #axis.text.x=element_blank(),
                              #axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                pQUESUB+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0.5,0.4,0.5,0),"cm")),
                nrow = 1,ncol = 3,label_size = 16,
                labels = c("a) FAGSYL","b) PINPINA","c) QUESUB"),label_x = 0,label_y = 1.025)
p1
p3 <- plot_grid(p1,legend,nrow=2,ncol=1,rel_heights = c(0.95,0.4),scale = c(0.95,0.05))

cairo_pdf(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/TestFigS14INTER.pdf"),width = 20,height = 1.75*4,family="helvetica")
p3
dev.off()
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/FigS12.ReviewV1.pdf"),plot = p3,base_width = 5.5, base_height = 1.75, dpi = 300,units = "in",nrow=4,ncol=3)




## S16
Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/S16/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
NAME <- c(1:22)

for (i in c(1:22)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}



p <- plot_grid(p1+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              # axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0,0,0,0),"cm")),
                
                p2+theme(legend.position = "none",
                               #axis.text.x=element_blank(),
                               #axis.text.y=element_blank(),
                               #axis.title.x=element_blank(),
                               axis.title.y=element_blank(),
                               #plot.title = element_blank(),
                               plot.margin = unit(c(0,0,0,0),"cm")),
                p3+theme(legend.position = "none",
                              # axis.text.x=element_blank(),
                              # axis.text.y=element_blank(),
                              #axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              #plot.title = element_blank(),
                              plot.margin = unit(c(0,0,0,0),"cm")),
                p4+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p5+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                
                p6+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p7+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p8+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p9+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                
                p10+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p11+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p12+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p13+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                
                p14+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p15+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p16+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p17+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                
                p18+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p19+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p20+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                p21+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                
                p22+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm")),
                legend,
                nrow = 6,ncol = 4,label_size = 16,
                labels = c("a) ACEPSE","b) ACEPSE","c) ALNGLU",
                           "d) BETPEN","e) PICABI","f) PICABI",
                           "g) PICABI","h) PICABI","i) PINNIG",
                           "j) PINNIG","k) PINPINA","l) PINPINA",
                           "m) PINPIN","n) PINPIN","o) PINSYL",
                           "p) PINSYL","q) POPNIG","r) QUEILE",
                           "s) QUEPYR","t) QUEROB","u) QUESUB","v) QUESUB"),rel_heights = c(0.95),scale = c(0.95),label_x = -0.05,label_y = 1)

p3
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/FigS16TEST.pdf"),plot = p,base_width = 8, base_height = 6, dpi = 300,units = "in",nrow=6,ncol=4)

## S17
Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/S17/") # Directory 
setwd(Dir)
#assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
NAME <- c(1:20)

for (i in c(1:20)){
  assign(paste0("p",NAME[i]),readRDS(list.files(pattern=".rds")[i]))
}


p <- plot_grid(p1+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        # axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               
               p2+theme(legend.position = "none",
                        #axis.text.x=element_blank(),
                        #axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p3+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p4+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p5+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        # axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               
               p6+theme(legend.position = "none",
                        #axis.text.x=element_blank(),
                        #axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p7+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p8+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p9+theme(legend.position = "none",
                        # axis.text.x=element_blank(),
                        # axis.text.y=element_blank(),
                        #axis.title.x=element_blank(),
                        # axis.title.y=element_blank(),
                        #plot.title = element_blank(),
                        plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               
               p10+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p11+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p12+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p13+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               
               p14+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p15+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p16+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p17+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         # axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               
               p18+theme(legend.position = "none",
                         #axis.text.x=element_blank(),
                         #axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p19+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               p20+theme(legend.position = "none",
                         # axis.text.x=element_blank(),
                         # axis.text.y=element_blank(),
                         #axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         #plot.title = element_blank(),
                         plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm")),
               nrow = 5,ncol = 4,label_size = 16,
               labels = c("a) ABIALB","b) ACEPSE","c) ACEPSE",
                          "d) ALNGLU","e) ALNGLU","f) BETPEN",
                          "g) FAGSYL","h) PICABI","i) PINNIG",
                          "j) PINNIG","k) PINPINA","l) PINPINA",
                          "m) PINSYL","n) PINSYL","o) POPTRE",
                          "p) QUEILE","q) QUEPYR","r) QUEPYR",
                          "s) QUEROB","t) QUROB"),label_x = -0.05,label_y = 1)

pp <- plot_grid(p,legend,nrow=2,ncol=1,rel_heights = c(1,0.05),scale = c(0.98,0.05))

cairo_pdf(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/TestFigS17INTER.pdf"),width = 8*4,height = 6*5,family="helvetica")
pp
dev.off()
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/FigS17TEST.pdf"),plot = p,base_width = 8, base_height = 6, dpi = 300,units = "in",nrow=5,ncol=4)







