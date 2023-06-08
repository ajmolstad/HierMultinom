setwd("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/Results/")
library(ggplot2)
library(reshape)
getMSE_Results <- function(paramSeq){
  Results <- matrix(0, nrow=4, ncol=length(paramSeq)*100)
  k <- 1
  rho <- NULL
  nreps <- 100
  for(kk in 1:length(paramSeq)){
    for(jj in 1:nreps){
      temp <- try(readRDS(paste(paramSeq[kk],"_", jj, ".RDS", sep="")))
      if(class(temp)!="try-error"){
        Results[,k] <- c(
          temp$glmnetGroup.KL,temp$glmnetL1.KL, temp$twoStepGrouped.KL, temp$Ours.KL)
        rho[k] <- strsplit(paramSeq[kk],"_p")[[1]][2]
        k <- k + 1
        cat(k, "\n")
        } 
    }
  }
  return(list(Results, rho))
}
greyCols <- c("grey100", "grey80", "grey60", "grey40")


paramSeq <- paste("Model1_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)



dat1 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat1$Method <- factor(dat1$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))


p1 <- ggplot(dat1, aes(x=as.factor(p), y=log(KL), fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 1") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model2_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)



dat2 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
                            t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
                   "Method" = c(rep("Group", length(t0[[2]])), 
                                rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
                   "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat2$Method <- factor(dat2$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))


p2 <- ggplot(dat2, aes(x=as.factor(p), y=log(KL), fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 2") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)")  + scale_fill_manual(values = greyCols) 






paramSeq <- paste("Model3_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat3 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
                            t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
                   "Method" = c(rep("Group", length(t0[[2]])), 
                                rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
                   "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat3$Method <- factor(dat3$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))

library(ggplot2)

p3 <- ggplot(dat3, aes(x=as.factor(p), y=log(KL), fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 3") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)")  + scale_fill_manual(values = greyCols) 




paramSeq <- paste("Model4_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat4 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
                            t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
                   "Method" = c(rep("Group", length(t0[[2]])), 
                                rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
                   "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat4$Method <- factor(dat4$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))


p4 <- ggplot(dat4, aes(x=as.factor(p), y=log(KL), fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 4") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model5_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat5 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
                            t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
                   "Method" = c(rep("Group", length(t0[[2]])), 
                                rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
                   "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat5$Method <- factor(dat5$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))


p5 <- ggplot(dat5, aes(x=as.factor(p), y=log(KL), fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 5") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model6_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat6 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
                            t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
                   "Method" = c(rep("Group", length(t0[[2]])), 
                                rep("L1", length(t0[[2]])), rep("Approx", length(t0[[2]])), rep("mrMLR", length(t0[[2]]))),
                   "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)
dat6$Method <- factor(dat6$Method, ordered=TRUE, levels = c("Group", "L1", "mrMLR", "Approx"))

p6 <- ggplot(dat6, aes(x=as.factor(p), y=log(KL), fill = Method)) + theme_bw() + 
  scale_fill_manual(values = greyCols) + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 6") + 
  theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)") 





library(cowplot)
pdf("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/KL_NOL.pdf", width=7.2, height=3.8)
plot_grid(p1 + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .05), "cm")), 
  p2 + ylab("") + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .0), "cm")), 
    p3 + ylab("") + xlab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.05, .05, .0, .0), "cm")), 
    p4 + ylab("log(KL divergence)") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .05), "cm")), 
    p5 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .0), "cm")),
    p6 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .05, .05, .0), "cm")), nrow=2, align='vh')
dev.off()

pdf("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/KL_NOL_Legend.pdf", width=6, height=3.6)
plot_grid(p6 + theme(legend.position="bottom"),nrow=1)
dev.off()












setwd("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/Results/")
library(ggplot2)
library(reshape)

getMSE_Results <- function(paramSeq){
  Results <- matrix(0, nrow=4, ncol=length(paramSeq)*100)
  k <- 1
  rho <- NULL
  nreps <- 100
  for(kk in 1:length(paramSeq)){
    for(jj in 1:nreps){
      temp <- try(readRDS(paste(paramSeq[kk],"_", jj, ".RDS", sep="")))
      if(class(temp)!="try-error"){
        Results[,k] <- c(
          temp$glmnetGroup.classErr,temp$glmnetL1.classErr, temp$twoStepGrouped.classErr, temp$Ours.classErr)
        rho[k] <- strsplit(paramSeq[kk],"_p")[[1]][2]
        k <- k + 1
        cat(k, "\n")
        } 
    }
  }
  return(list(Results, rho))
}

paramSeq <- paste("Model1_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat1 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

  library(ggplot2)

p1 <- ggplot(dat1, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 1") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Classification error") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model2_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)



dat2 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)

p2 <- ggplot(dat2, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 2") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Classification error") + scale_fill_manual(values = greyCols) 






paramSeq <- paste("Model3_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat3 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)

p3 <- ggplot(dat3, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 3") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Classification error") + scale_fill_manual(values = greyCols) 




paramSeq <- paste("Model4_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat4 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p4 <- ggplot(dat4, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 4") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Classification error") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model5_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat5 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p5 <- ggplot(dat5, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 5") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("log(KL divergence)") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model6_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat6 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p6 <- ggplot(dat6, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 6") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Classification error") + scale_fill_manual(values = greyCols) 




library(cowplot)
pdf("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/CE_NOL.pdf", width=7.2, height=3.8)
plot_grid(p1 + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .05), "cm")), 
  p2 + ylab("") + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .0), "cm")), 
    p3 + ylab("") + xlab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.05, .05, .0, .0), "cm")), 
    p4 + ylab("Classification error") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .05), "cm")), 
    p5 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .0), "cm")),
    p6 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .05, .05, .0), "cm")), nrow=2, align='vh')
dev.off()




















setwd("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/Results/")
library(ggplot2)
library(reshape)

getMSE_Results <- function(paramSeq){
  Results <- matrix(0, nrow=4, ncol=length(paramSeq)*100)
  k <- 1
  rho <- NULL
  nreps <- 100
  for(kk in 1:length(paramSeq)){
    for(jj in 1:nreps){
      temp <- try(readRDS(paste(paramSeq[kk],"_", jj, ".RDS", sep="")))
      if(class(temp)!="try-error"){
        Results[,k] <- c(
          temp$glmnetGroup.hellinger,temp$glmnetL1.hellinger, temp$twoStepGrouped.hellinger, temp$Ours.hellinger)
        rho[k] <- strsplit(paramSeq[kk],"_p")[[1]][2]
        k <- k + 1
        cat(k, "\n")
        } 
    }
  }
  return(list(Results, rho))
}

paramSeq <- paste("Model1_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat1 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

  library(ggplot2)

p1 <- ggplot(dat1, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 1") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model2_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)



dat2 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)

p2 <- ggplot(dat2, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 2") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 






paramSeq <- paste("Model3_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat3 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))

library(ggplot2)

p3 <- ggplot(dat3, aes(x=as.factor(p), y=KL, fill = as.factor(Method))) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 3") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 




paramSeq <- paste("Model4_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat4 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p4 <- ggplot(dat4, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 4") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model5_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat5 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p5 <- ggplot(dat5, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 5") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 



paramSeq <- paste("Model6_p", c(100, 200, 500, 1000), sep="")
t0 <- getMSE_Results(paramSeq)


dat6 <- data.frame("KL" = c(t0[[1]][1,1:length(t0[[2]])], 
    t0[[1]][2,1:length(t0[[2]])], t0[[1]][3,1:length(t0[[2]])], t0[[1]][4,1:length(t0[[2]])]), 
    "Method" = c(rep("Group", length(t0[[2]])), 
      rep("L1", length(t0[[2]])), rep("Two-Step", length(t0[[2]])), rep("mrMLRolution", length(t0[[2]]))),
    "p" = rep(as.numeric(t0[[2]]), 4))


library(ggplot2)

p6 <- ggplot(dat6, aes(x=as.factor(p), y=KL, fill = Method)) + theme_bw() + 
  geom_boxplot(lwd=.5,fatten=.65, position=position_dodge(.65),  outlier.size = .4) + xlab("p") + ggtitle("Model 6") + theme(plot.title = element_text(hjust = 0.5, size=10)) + 
  ylab("Hellinger distance") + scale_fill_manual(values = greyCols) 




library(cowplot)
pdf("/blue/amolstad/amolstad/HierMultinom/Simulations_R1/Hellinger_NOL.pdf", width=7.2, height=3.8)
plot_grid(p1 + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .05), "cm")), 
  p2 + ylab("") + xlab("") + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .0), "cm")), 
    p3 + ylab("") + xlab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.05, .05, .0, .0), "cm")), 
    p4 + ylab("Hellinger distance") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .05), "cm")), 
    p5 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .0), "cm")),
    p6 + ylab("") + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .05, .05, .0), "cm")), nrow=2, align='vh')
dev.off()



