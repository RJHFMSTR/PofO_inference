setwd('~/Dropbox/Parent-of-Origin/Manuscript/15_09_2022-NC_revision_2/Figures/Source_Data/')
library(ggplot2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggtext)
library(glue)
mycols<-brewer.pal(9, 'Set1')

a4_width=210 # mm
a4_height=297 # mm
line_size=0.8
point_size=1.5
st=8 # size lab text
stt=5
lt=6
sig<-5e-08


# Figure 5A
  
  d<-read.table('source_data_F5A.txt', hea=T)
  
  g_ss<-ggplot(d, aes(x=N, y=median_log10_p, group=mode, color=mode)) + 
    geom_line() +
    geom_point() + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
    geom_errorbar(aes(x=N,ymin=median_log10_p-std_log10_p, ymax=median_log10_p+std_log10_p), width=5,
                  position=position_dodge(0.2)) + theme_bw() +
    theme(legend.position = c(0.006,.994),legend.justification = c(0, 1),legend.direction="vertical",
          legend.box.background = element_rect(color="black", size=1),
          panel.border = element_rect(color="black", size=1),
          axis.ticks = element_line(size=.7), axis.text = element_text(size=st),
          axis.title = element_text(size=st), legend.text = element_text(size=lt),
          legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
          legend.key.width = unit(.1,'cm'), legend.key.size = unit(0.5, 'cm'),
          legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.05,'cm'),
          title = element_text(size=stt), plot.title = element_text(face='bold',hjust = 0.5),
          
          panel.grid = element_blank())+ 
    guides(col=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=.4))) +
    labs(y='-log10(p-values)', x='Number of samples', col='P-values :') +
    scale_x_continuous(breaks=c(5000,10000,15000,20000,26393, 30000),labels=c(4909,10000,15000,20000,26393, 30000)) +
    scale_color_manual(values=c(mycols[9], mycols[3], mycols[1], mycols[2]),
                       labels=c('additive', 'maternal','paternal','differential')) + lims(y=c(0,ylimm))

  
  title <- ggdraw() + 
    draw_label('',
               fontface = 'bold', hjust = 0.45  , size=st) +  theme(plot.margin = margin(0, 0, -5, 7)  ) +
    theme(plot.background = element_rect(fill='white', colour = 'white'))
  
  
  Figure_5A<-plot_grid(title,g_ss, ncol=1,rel_heights = c(0.08, 1))
  

  # Figure 5B
  
  d<-read.table('source_data_F5B.txt', hea=T)
  
  g_ss<-ggplot(d, aes(x=N, y=median_log10_p, group=mode, color=mode)) + 
    geom_line() +
    geom_point() + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
    geom_errorbar(aes(x=N,ymin=median_log10_p-std_log10_p, ymax=median_log10_p+std_log10_p), width=5,
                  position=position_dodge(0.2)) + theme_bw() +
    theme(legend.position = c(0.006,.994),legend.justification = c(0, 1),legend.direction="vertical",
          legend.box.background = element_rect(color="black", size=1),
          panel.border = element_rect(color="black", size=1),
          axis.ticks = element_line(size=.7), axis.text = element_text(size=st),
          axis.title = element_text(size=st), legend.text = element_text(size=lt),
          legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
          legend.key.width = unit(.1,'cm'), legend.key.size = unit(0.5, 'cm'),
          legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.05,'cm'),
          title = element_text(size=stt), plot.title = element_text(face='bold',hjust = 0.5),
          
          panel.grid = element_blank())+ 
    guides(col=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=.4))) +
    labs(y='-log10(p-values)', x='Number of samples', col='P-values :') +
    scale_x_continuous(breaks=c(5000,10000,15000,20000,26393, 30000),labels=c(4909,10000,15000,20000,26393, 30000)) +
    scale_color_manual(values=c(mycols[9], mycols[3], mycols[1], mycols[2]),
                       labels=c('additive', 'maternal','paternal','differential')) + lims(y=c(0,ylimm))
  
  
  title <- ggdraw() + 
    draw_label('',
               fontface = 'bold', hjust = 0.45  , size=st) +  theme(plot.margin = margin(0, 0, -5, 7)  ) +
    theme(plot.background = element_rect(fill='white', colour = 'white'))
  
  
  Figure_5B<-plot_grid(title,g_ss, ncol=1,rel_heights = c(0.08, 1))
  

# Figure 5C
  d<-as.data.frame(data.table::fread('source_data_F5C.txt', hea=T))
  ylimm<-max(d$median_log10_p)+1
  g_r<-ggplot(d, aes(x=N, y=median_log10_p, group=mode, color=mode)) + 
    geom_line() +
    geom_point() +
    geom_errorbar(aes(x=N,ymin=median_log10_p-std_log10_p, ymax=median_log10_p+std_log10_p), width=5,
                  position=position_dodge(0.2)) + theme_bw() +
    geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
    theme(legend.position = c(.992,.993),legend.justification = c(1, 1),legend.direction="vertical",
          legend.box.background = element_rect(color="black", size=1),
          panel.border = element_rect(color="black", size=1),
          axis.ticks = element_line(size=.7), axis.text = element_text(size=st),
          axis.title = element_text(size=st), legend.text = element_text(size=lt),
          legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
          legend.key.width = unit(.1,'cm'), legend.key.size = unit(0.5, 'cm'),
          legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.01,'cm'),
          title = element_text(size=stt), plot.title = element_text(face='bold',hjust = 0.5),
          
          panel.grid = element_blank())+ 
    guides(col=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=.4))) +
    labs(y='-log10(p-values)', x='% of PofO errors', col='P-values :') +
    scale_x_continuous(breaks=c(100,400,1000,2000,6000,12000), labels=c(0,2,5,10,25,50)) +
    scale_color_manual(values=c(mycols[9], mycols[3], mycols[1], mycols[2]),
                       labels=c('additive', 'maternal','paternal','differential')) + lims(y=c(0,ylimm))
  g_r
  title <- ggdraw() + 
    draw_label('',
               fontface = 'bold', hjust = 0.45  , size=st) +  theme(plot.margin = margin(0, 0, -5, 7)  ) +
    theme(plot.background = element_rect(fill='white', colour = 'white'))
  
  
  Figure_5C<-plot_grid(title,g_r, ncol=1,rel_heights = c(0.08, 1))
  
  
  
  # Figure 5D
  d<-as.data.frame(data.table::fread('source_data_F5D.txt', hea=T))
  ylimm<-max(d$median_log10_p)+1
  g_r<-ggplot(d, aes(x=N, y=median_log10_p, group=mode, color=mode)) + 
    geom_line() +
    geom_point() +
    geom_errorbar(aes(x=N,ymin=median_log10_p-std_log10_p, ymax=median_log10_p+std_log10_p), width=5,
                  position=position_dodge(0.2)) + theme_bw() +
    geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
    theme(legend.position = c(.992,.993),legend.justification = c(1, 1),legend.direction="vertical",
          legend.box.background = element_rect(color="black", size=1),
          panel.border = element_rect(color="black", size=1),
          axis.ticks = element_line(size=.7), axis.text = element_text(size=st),
          axis.title = element_text(size=st), legend.text = element_text(size=lt),
          legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
          legend.key.width = unit(.1,'cm'), legend.key.size = unit(0.5, 'cm'),
          legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.01,'cm'),
          title = element_text(size=stt), plot.title = element_text(face='bold',hjust = 0.5),
          
          panel.grid = element_blank())+ 
    guides(col=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=.4))) +
    labs(y='-log10(p-values)', x='% of PofO errors', col='P-values :') +
    scale_x_continuous(breaks=c(100,400,1000,2000,6000,12000), labels=c(0,2,5,10,25,50)) +
    scale_color_manual(values=c(mycols[9], mycols[3], mycols[1], mycols[2]),
                       labels=c('additive', 'maternal','paternal','differential')) + lims(y=c(0,ylimm))
  title <- ggdraw() + 
    draw_label('',
               fontface = 'bold', hjust = 0.45  , size=st) +  theme(plot.margin = margin(0, 0, -5, 7)  ) +
    theme(plot.background = element_rect(fill='white', colour = 'white'))
  
  
  Figure_5D<-plot_grid(title,g_r, ncol=1,rel_heights = c(0.08, 1))
  
  
Figure5<-plot_grid(Figure_5A, Figure_5B, Figure_5C,Figure_5D,ncol=2)


#ggsave(file="Final_figures/Figure5.png", assembly, width = a4_width, height = a4_height/1.5, units = "mm")
