#================================================================================================================
# Plots for the simulation study
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 10-January-2024
#================================================================================================================

library(dplyr)
library(latex2exp)
library(ggplot2)
library(ggpubr) 

# figures for b_p
# to load the results from "df_bp_R1" data frame
load("df_bp_R1.RData")

labels.CorrXpZp <- c("0" = "Correlation 0", "0.2" = "Correlation 0.2", "0.5" = "Correlation 0.5")
labels.SNR <- c("0.5" = "v = 0.5", "1" = "v = 1", "1.5"= "v = 1.5", "2"= "v = 2","2.5"= "v = 2.5" )
labels.bp.methods = expression( paste("penMCFM(",alpha, "=0.1)"), 
                                paste("penMCFM(",alpha, "=0.5)"),
                                paste("penMCFM(",alpha, "=0.9)"),
                                paste("penMCFM(",alpha, "=1)"),
                                "penMCFM(GMIFS)","MCM(GMIFS)")

# Line types and colors for c("penMCFM (EM) 0.1", "penMCFM (EM) 0.5","penMCFM (EM) 0.9","penMCFM (EM) 1","penMCFM (GMIFS)", "MCM (GMIFS)") ), respectively.
linetype.bp = c("solid","solid","solid","dotted","solid", "solid", "dotted") 
color.bp = c('skyblue4','darkorange','coral3','black','darkolivegreen4','azure4' ) 

df.bp.nonzeros.R1 = df_bp_R1 %>% 
  group_by(SNR, CorrXpZp, Method) %>% 
  summarise(y=mean(nonzeros_bp), sd=sd(nonzeros_bp), median=median(nonzeros_bp), n=n() )

bp.nonzeros.R1 <- ggplot(df.bp.nonzeros.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.bp, labels = labels.bp.methods )+
  scale_color_manual(values = color.bp,  labels = labels.bp.methods )+
  ylab(TeX(r'(Number of nonzeros for $b_{p}$)'))+
  xlab(TeX(r'($v$)') )  +
  geom_hline(yintercept = 20, linetype="dotted")+
  guides( color=guide_legend(nrow=1, ncol=6, title=NULL), group=guide_legend(nrow=1, ncol=6, title=NULL),
          linetype=guide_legend(nrow=1, ncol=6, title=NULL))

plot.bp.nonzeros.R1 = bp.nonzeros.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller(CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7))
plot.bp.nonzeros.R1


# sensitivity of b_p
df.bp.Sensitivity.R1 = df_bp_R1 %>% 
  group_by(SNR, CorrXpZp, Method) %>% 
  summarise(y=mean(Sensitivity_bp), sd=sd(Sensitivity_bp), media=median(Sensitivity_bp), n=n() )

bp.Sensitivity.R1 <- ggplot(df.bp.Sensitivity.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position = position_dodge(0.5), linewidth = 1) +
  geom_point(position = position_dodge(0.5))+
  scale_linetype_manual(values = linetype.bp, labels = labels.bp.methods )+
  scale_color_manual(values = color.bp, labels = labels.bp.methods )+
  ylab(TeX( r'(Sensitivity for $b_{p}$)'))+
  xlab(TeX(r'($v$)') )  +
  guides( color = guide_legend(nrow=1, ncol=6, title=NULL), group = guide_legend(nrow=1, ncol=6, title=NULL),
          linetype = guide_legend(nrow=1, ncol=6, title=NULL))

plot.bp.Sensitivity.R1 = bp.Sensitivity.R1 + 
  facet_grid(~ CorrXpZp, labeller = labeller(CorrXpZp=labels.CorrXpZp) )+
  theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) )
plot.bp.Sensitivity.R1


# Specificity of b_p
df.bp.Specificity.R1 = df_bp_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(Specificity_bp), sd=sd(Specificity_bp), median=median(Specificity_bp), n=n() )

bp.Specificity.R1 <- ggplot(df.bp.Specificity.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.bp, labels = labels.bp.methods )+
  scale_color_manual(values = color.bp, labels = labels.bp.methods )+
  ylab(TeX( r'(Specificity for $b_{p}$)')) +
  xlab(TeX(r'($v$)') ) +
  guides( color = guide_legend(nrow=1, ncol=6, title=NULL), group = guide_legend(nrow=1, ncol=6, title=NULL),
          linetype = guide_legend(nrow=1, ncol=6, title=NULL))

plot.bp.Specificity.R1 = bp.Specificity.R1 + 
  facet_grid( ~ CorrXpZp, labeller = labeller( CorrXpZp = labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7))
plot.bp.Specificity.R1


# FPR of b_p
df.bp.FPR.R1 = df_bp_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(FPR_bp), sd=sd(FPR_bp), median=median(FPR_bp), n=n() )

bp.FPR.R1 <- ggplot(df.bp.FPR.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.bp, labels = labels.bp.methods )+
  scale_color_manual(values = color.bp, labels = labels.bp.methods )+
  ylab(TeX( r'(FPR for $b_{p}$)')) + 
  xlab(TeX(r'($v$)') ) +
  guides( color = guide_legend(nrow=1, ncol=6, title=NULL), group = guide_legend(nrow=1, ncol=6, title=NULL),
          linetype = guide_legend(nrow=1, ncol=6, title=NULL))

plot.bp.FPR.R1 = bp.FPR.R1 + 
  facet_grid(~ CorrXpZp, labeller = labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() + theme(legend.position="bottom", legend.text = element_text(size = 7) )
plot.bp.FPR.R1 


gg.bp.selection.performance = ggarrange(plot.bp.nonzeros.R1, plot.bp.Sensitivity.R1, plot.bp.Specificity.R1, plot.bp.FPR.R1,
                                        common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 4 ) 
gg.bp.selection.performance

ggsave(file = "bp.selection.performance.plot.pdf", plot = gg.bp.selection.performance, width = 210, height = 297, units = "mm", dpi= 360 )


# Absolute bias of uncured estimation 
df.uncure.bias.R1 = df_bp_R1 %>%  group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(uncure.bias), sd=sd(uncure.bias), median=median(uncure.bias), n=n() )

bp.uncure.bias.R1 <- ggplot(df.uncure.bias.R1, aes(x=SNR, y=abs(y), linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.bp, labels = labels.bp.methods)+
  scale_color_manual(values = color.bp, labels = labels.bp.methods)+
  ylab(TeX( r'(Absolute Bias of $\widehat{\pi}(z)$)' ) ) +
  xlab(TeX(r'($v$)') ) + 
  guides( color=guide_legend(nrow=1, ncol=6, title=NULL), group=guide_legend(nrow=1, ncol=6, title=NULL),
          linetype=guide_legend(nrow=1, ncol=6, title=NULL)) 

plot.uncure.bias.R1 = bp.uncure.bias.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom",  legend.text = element_text(size = 7) )
plot.uncure.bias.R1

# MSE of uncured estimation
df.uncure.mse.R1 = df_bp_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(uncure.mse), sd=sd(uncure.mse), median=median(uncure.mse), n=n() )

bp.uncure.mse.R1 <- ggplot(df.uncure.mse.R1, aes(x=SNR, y=y, linetype=Method,  group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values=linetype.bp, labels = labels.bp.methods)+
  scale_color_manual(values=color.bp, labels = labels.bp.methods)+
  ylab(TeX( r'(MSE of $\widehat{\pi}(z)$)' ) ) + 
  xlab(TeX(r'($v$)') )+ 
  guides( color=guide_legend(nrow=1, ncol=6, title=NULL), group=guide_legend(nrow=1, ncol=6, title=NULL),
          linetype=guide_legend(nrow=1, ncol=6, title=NULL)) 

plot.uncure.mse.R1 = bp.uncure.mse.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom",  legend.text = element_text(size = 7) )
plot.uncure.mse.R1


gg.uncured.rate.bias.mse = ggarrange(plot.uncure.bias.R1, plot.uncure.mse.R1, 
                                     common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 2 ) 
gg.uncured.rate.bias.mse

ggsave(file = "uncured.rate.bias.mse.plot.pdf", plot = gg.uncured.rate.bias.mse, width = 210, height = 150, units = "mm", dpi= 360 )


# figures for beta_p
# to load the results from "df_betap_R1" data frame
load("df_betap_R1.RData")

labels.betap.methods = expression( paste("penMCFM(",alpha, "=0.1)"), 
                                   paste("penMCFM(",alpha, "=0.5)"),
                                   paste("penMCFM(",alpha, "=0.9)"),
                                   paste("penMCFM(",alpha, "=1)"),
                                   "penMCFM(GMIFS)","MCM(GMIFS)", "penCox.1se")

# Line types and colors for c("penMCFM (EM) 0.1", "penMCFM (EM) 0.5","penMCFM (EM) 0.9","penMCFM (EM) 1","penMCFM (GMIFS)", "MCM (GMIFS)", "penCox.1se") ), respectively.
linetype.betap = c("solid","solid","solid","dotted","solid", "solid", "dotted") 
color.betap = c('skyblue4','darkorange','coral3','black','darkolivegreen4','azure4','magenta3') 


# non-zeros beta_p
df.betap.nonzeros.R1 = df_betap_R1 %>% 
  group_by(SNR, CorrXpZp, Method) %>% 
  summarise(y=mean(nonzeros_betap), sd=sd(nonzeros_betap), median=median(nonzeros_betap), n=n() )

betap.nonzeros.R1 <- ggplot(df.betap.nonzeros.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.betap, labels = labels.betap.methods )+
  scale_color_manual(values = color.betap, labels = labels.betap.methods )+
  ylab(TeX(r'(Number of nonzeros for $\beta_{p}$)'))+
  xlab(TeX(r'($v$)') )  +
  geom_hline(yintercept = 20, linetype="dotted")+
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL ), group=guide_legend(nrow=1, ncol=7, title=NULL ),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL) )

plot.betap.nonzeros.R1 = betap.nonzeros.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller(CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7) )
plot.betap.nonzeros.R1


# sensitivity of beta_p
df.betap.Sensitivity.R1 = df_betap_R1 %>% 
  group_by(SNR, CorrXpZp, Method) %>% 
  summarise(y=mean(Sensitivity_betap), sd=sd(Sensitivity_betap), media=median(Sensitivity_betap), n=n() )

betap.Sensitivity.R1 <- ggplot(df.betap.Sensitivity.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.betap, labels = labels.betap.methods )+
  scale_color_manual(values = color.betap, labels = labels.betap.methods )+
  ylab(TeX( r'(Sensitivity for $\beta_{p}$)'))+
  xlab(TeX(r'($v$)') )  +
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL), group=guide_legend(nrow=1, ncol=7, title=NULL),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL))

plot.betap.Sensitivity.R1 = betap.Sensitivity.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller(CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7) ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) )
plot.betap.Sensitivity.R1


# Specificity of beta_p
df.betap.Specificity.R1 = df_betap_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(Specificity_betap), sd=sd(Specificity_betap), median=median(Specificity_betap), n=n() )

betap.Specificity.R1 <- ggplot(df.betap.Specificity.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.betap, labels = labels.betap.methods )+
  scale_color_manual(values = color.betap, labels = labels.betap.methods )+
  ylab(TeX( r'(Specificity for $\beta_{p}$)')) +
  xlab(TeX(r'($v$)') ) +
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL), group=guide_legend(nrow=1, ncol=7, title=NULL),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL))

plot.betap.Specificity.R1 = betap.Specificity.R1 + 
  facet_grid( ~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7) )
plot.betap.Specificity.R1


## FPR plot for beta_p
df.betap.FPR.R1 = df_betap_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(FPR_betap), sd=sd(FPR_betap), median=median(FPR_betap), n=n() )

betap.FPR.R1 <- ggplot(df.betap.FPR.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.betap, labels = labels.betap.methods )+
  scale_color_manual(values = color.betap, labels = labels.betap.methods )+
  ylab(TeX( r'(FPR for $\beta_{p}$)')) + 
  xlab(TeX(r'($v$)') ) +
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL), group=guide_legend(nrow=1, ncol=7, title=NULL),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL) )

plot.betap.FPR.R1 = betap.FPR.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() + theme(legend.position="bottom", legend.text = element_text(size = 7) )  
plot.betap.FPR.R1 

gg.betap.selection.performance = ggarrange(plot.betap.nonzeros.R1, plot.betap.Sensitivity.R1, plot.betap.Specificity.R1, plot.betap.FPR.R1,
                                           common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 4 ) 
gg.betap.selection.performance

ggsave(file = "betap.selection.performance.plot.pdf", plot = gg.betap.selection.performance, width = 210, height = 297, units = "mm", dpi= 360 )


# C-statistics plot based on train data from the data frame "df_bp_R1"
df.Cstat.train.R1 = df_bp_R1 %>% group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(Cstat), sd=sd(Cstat), median=median(Cstat), n=n() )

# to load the penCox.1se.train results from "df.penCox.1se.train.Cstat" data frame
load("df.penCox.1se.train.Cstat.RData")

df.Cstat.train.R1 = rbind(df.Cstat.train.R1, df.penCox.1se.train.Cstat)
df.Cstat.train.R1$Method = as.factor( df.Cstat.train.R1$Method )
labels.bp.methods.R1 = expression( paste("penMCFM(",alpha, "=0.1)"), 
                                   paste("penMCFM(",alpha, "=0.5)"),
                                   paste("penMCFM(",alpha, "=0.9)"),
                                   paste("penMCFM(",alpha, "=1)"),
                                   "penMCFM(GMIFS)","MCM(GMIFS)", "penCox.1se")

bp.Cstat.train.R1 <- ggplot(df.Cstat.train.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values=c("solid","solid","solid","dotted","solid", "solid", "dotted"),  labels = labels.bp.methods.R1)+
  scale_color_manual(values=c('skyblue4','darkorange','coral3','black','darkolivegreen4','azure4','magenta3'),  labels = labels.bp.methods.R1)+
  ylab("C-statistic") +
  xlab(TeX(r'($v$)') )+ 
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL), group=guide_legend(nrow=1, ncol=7, title=NULL),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL)) 

plot.Cstat.train.R1 = bp.Cstat.train.R1 + 
  facet_grid(~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) )+
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7) ) 
plot.Cstat.train.R1


# Cstat test set plot from the data frame "df_betap_R1"
df.Cstat.test.R1 = df_betap_R1 %>%
  group_by(SNR, CorrXpZp, Method) %>%
  summarise(y=mean(Cstat_test), sd=sd(Cstat_test), median=median(Cstat_test), n=n() )

Cstat.test.R1 <- ggplot(df.Cstat.test.R1, aes(x=SNR, y=y, linetype=Method, group=Method, color=Method)) + 
  geom_line(position=position_dodge(0.5), linewidth=1) +
  geom_point(position=position_dodge(0.5))+
  scale_linetype_manual(values = linetype.betap, labels = labels.betap.methods )+
  scale_color_manual(values = color.betap, labels = labels.betap.methods )+
  ylab("C-statistic") + 
  xlab( TeX(r'($v$)') ) +
  guides( color=guide_legend(nrow=1, ncol=7, title=NULL), group=guide_legend(nrow=1, ncol=7, title=NULL),
          linetype=guide_legend(nrow=1, ncol=7, title=NULL) )

plot.Cstat.test.R1 = Cstat.test.R1 + 
  facet_grid( ~ CorrXpZp, labeller=labeller( CorrXpZp=labels.CorrXpZp) ) +
  theme_bw() +  theme(legend.position="bottom", legend.text = element_text(size = 7) )
plot.Cstat.test.R1

gg.Cstat.train.test = ggarrange(plot.Cstat.train.R1, plot.Cstat.test.R1,
                                common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 2 ) 
gg.Cstat.train.test

ggsave(file = "Cstat.train.test.plot.pdf", plot = gg.Cstat.train.test, width = 210, height = 150, units = "mm", dpi= 360 )
