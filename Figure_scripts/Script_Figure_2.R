library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
mycols<-brewer.pal(9, 'Set1')

a4_width=210 # mm
a4_height=297 # mm
line_size=0.8
point_size=1.5
st=6 # size lab text
stt=5
title_size=5
lt=5


#Figure 2A.
d<-as.data.frame(data.table::fread('source_data_F2A.txt', hea=T))
gA<-ggplot(d, aes(x=call_rate, y=error_rate, col=factor(cm), group=factor(cm))) + geom_point(size=point_size-0.8) + geom_line(size=line_size-0.3) +
  scale_color_manual(values=mycols[1:6]) + theme_bw() + theme(legend.position = c(0.006,.994),legend.justification = c(0, 1),legend.direction="horizontal",
                                                              legend.box.background = element_rect(color="black", size=1),
                                                              panel.border = element_rect(color="black", size=1),
                                                              axis.ticks = element_line(size=.4), axis.text = element_text(size=st),
                                                              axis.title = element_text(size=st), legend.text = element_text(size=lt),
                                                              legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
                                                              legend.key.width = unit(0,'cm'), legend.key.size = unit(0, 'cm'),
                                                              legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0,'cm'),
                                                              title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5),
                                                              plot.margin = margin(1, 5.5, 5.5, 9),
                                                              axis.title.y = element_text(vjust = 3.5),
                                                              
                                                              
                                                              panel.grid = element_blank())+ 
  guides(col=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size=.5))) +
  labs(x='Call rate [%]', y='Error rate [%]', col='Length of haplotype\n    segment [cM]') +
  geom_segment(aes(x = 85, y = 0.25, xend = d$call_rate[d$cm==3 & d$threshold==0.7]+1.5, yend = d$error_rate[d$cm==3 & d$threshold==0.7]-0.05), 
               arrow = arrow(length = unit(0.2, "cm")), col='darkgrey', size=0.2)


title <- ggdraw() + 
  draw_label('A. PofO inference parameters optimization',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))






Figure_2A<-plot_grid(title, gA, ncol=1,rel_heights = c(0.03,1))


#Figure 2B.

db<-as.data.frame(data.table::fread('source_data_F2B.txt', hea=T))
gB <- ggplot(db, aes(x = factor(pattern))) + geom_line(aes(y = error_rate, group=1, col='1'), size=line_size-0.3) +
  geom_line(aes(y = call_rate/100, group=1, col='2'), size=line_size-0.3)+geom_point(aes(y = error_rate, group=1, col='1'), size=point_size-0.8) +
  geom_point(aes(y = call_rate/100, group=1, col='2'), size=point_size-0.8)+
  scale_y_continuous(sec.axis = sec_axis(~.*100, name='Call rate [%]')) +
  labs(y='Error rate [%]', x='Group pattern') + scale_x_discrete(limits=c("0/1","0/2","1/1","1/2",">2"), labels=c(' ',' ',' ',' ',' ')) +
  theme_bw() + scale_color_manual(values=c(mycols[1], mycols[3]), label=c('Error rate','Call rate')) +  labs(x='', col='')  +
  theme(legend.position = c(0.008,.994),legend.justification = c(0, 1),legend.direction="horizontal",
        legend.box.background = element_rect(color="black", size=1),
        panel.border = element_rect(color="black", size=1),
        axis.ticks = element_line(size=.4), axis.text = element_text(size=st),
        axis.title = element_text(size=st), legend.text = element_text(size=lt),
        legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
        legend.key.width = unit(0,'cm'), legend.key.size = unit(0, 'cm'),
        legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.04,'cm'),
        title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5),
        plot.margin = margin(1, 5.5, 5.5, 7), panel.grid = element_blank()
        #axis.title.y.left = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        #axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
  )+ 
  guides(col=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size=.5))) 


title <- ggdraw() + 
  draw_label('B. PofO inference accuracy with optimal parameters',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))






Figure_2B<-plot_grid(title, gB, ncol=1,rel_heights = c(0.03,1))




# Figure 2C
pattern<-as.data.frame(data.table::fread('source_data_F2C.txt', hea=T))
gC<-ggplot(pattern, aes(fill=group, y=samples, x=factor(pattern, levels=c('0/1','0/2','1/1','1/2','>2')))) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c('grey','black'), labels=c('Validation cohort (N=1,399)','Full cohort (N=21,484)')) +
  labs(x='Group pattern', y='% of samples') +
  theme_bw()+  geom_text(aes(label=label), vjust=-0.3, color="black",
                         position = position_dodge(1), size=1.5) +
  theme(legend.position = c(.992,.993),legend.justification = c(1, 1),legend.direction="horizontal",
        legend.box.background = element_rect(color="black", size=1),
        panel.border = element_rect(color="black", size=1),
        axis.ticks = element_line(size=.4), axis.text = element_text(size=st),
        axis.title = element_text(size=st), legend.text = element_text(size=lt),
        legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
        legend.key.width = unit(.2,'cm'), legend.key.size = unit(1, 'cm'),
        legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0.05,'cm'),
        title = element_text(size=stt), plot.title = element_text(face='bold',hjust = 0.5),
        plot.margin = margin(1, 5.5, 5.5, 7),
        
        
        panel.grid = element_blank())+ 
  guides(fill=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size=.5))) +
  scale_x_discrete(limits=c("0/1","0/2","1/1","1/2",">2"), labels=c(' ',' ',' ',' ',' ')) + labs(x='', fill='') 


title <- ggdraw() + 
  draw_label('C. Surrogate parent group composition across cohorts',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))




Figure_2C<-plot_grid(title, gC, ncol=1,rel_heights = c(0.03,1))







title <- ggdraw() + 
  draw_label('Figure 2. Validation of the PofO inference',
             fontface = 'bold',x = 0,hjust = 0  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))


assembly_ABC<-plot_grid(Figure_2A, Figure_2B, Figure_2C, ncol=3, rel_widths=c(1,1.1,1))
assembly_ABC2<-plot_grid( title, assembly_ABC, nrow=2, rel_heights = c(0.05,1))



# Figure 2D

d<-as.data.frame(data.table::fread('source_data_F2DE.txt', hea=T))
d$error_rate[d$error_rate>10]<-11
gd1<-ggplot(d, aes(x=norm_pos, y=error_rate)) + geom_point(size=point_size/2,alpha=.1, shape=20) + geom_smooth(col='red') +
  theme_bw() +
  theme(panel.border = element_rect(color="black", size=1),
        axis.ticks = element_blank(), axis.text = element_text(size=st),
        axis.title.y  = element_text(size=st), legend.text = element_text(size=4),axis.title.x  = element_text(size=st-2),
        title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5), plot.margin=unit(c(1,5.5,-7.0,5.5), "pt"),
        panel.grid = element_blank())+ 
  guides(col=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size=.5))) +
  labs(x='', y='Error rate [%]', col='') + scale_x_continuous(limits=c(0,1), breaks=c(0,1), labels=c('','')) +
  geom_hline(yintercept=10.8, col='grey', linetype='dashed') + scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=c(0,2,4,6,8,10))


gd2<-ggplot(d, aes(x=norm_pos, y=call_rate)) + geom_point(size=point_size/2,alpha=.1, shape=20) + geom_smooth(col='red') +
  theme_bw() +
  theme(panel.border = element_rect(color="black", size=1),
        axis.ticks = element_line(size=.4), axis.text.x = element_text(size=st-2),axis.text.y = element_text(size=st),
        axis.title = element_text(size=st), legend.text = element_text(size=4),
        title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5), plot.margin=unit(c(-7.0,5.5,5.5,5.5), "pt"),
        
        
        panel.grid = element_blank())+ 
  guides(col=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size=.5))) +
  labs(x='Normalized chromosomes positions', y='Call rate [%]', col='')  + scale_x_continuous(limits=c(0,1), breaks=c(0,1),
                                                                                              labels=c('Start','End'))

title <- ggdraw() + 
  draw_label('D. Density of rates along chromosomes',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))

Figure_2D<-plot_grid(title, gd1,gd2, ncol=1,rel_heights = c(0.05, 1,1))



# 2E
d<-as.data.frame(data.table::fread('source_data_F2DE.txt', hea=T))
ge<-ggplot(d) + geom_histogram(aes(error_rate), binwidth = .4) + scale_y_log10(breaks=c(1,2,5,10,100, 1000, 10000, 100000),
                                                                               labels=c('1','2','5','10','100','1e3','1e4','1e5'),
                                                                               sec.axis = sec_axis(~.*100, name='', labels=rep('', 4)))+
  theme_bw() +
  theme(panel.border = element_rect(color="black", size=1),
        axis.ticks = element_line(size=.4), axis.text = element_text(size=st),
        axis.title = element_text(size=st), legend.text = element_text(size=4),
        title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5),
        axis.text.y=element_text(), axis.ticks.y.right  = element_blank(),
        plot.margin = margin(1, 12.5, 5.5, 5.5),
        
        
        panel.grid = element_blank())+ 
  guides(col=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size=.5)))+
  labs(x='Error rate [%]', y='# of variants (log10)', col='') 

title <- ggdraw() + 
  draw_label('E. Distribution of the number of variants per error rate',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))

Figure_2E<-plot_grid(title, ge, ncol=1,rel_heights = c(0.03, 1))



# 3F
d<-as.data.frame(data.table::fread('source_data_F2F.txt', hea=T))


# Annotation layer: add the chr label next to the chromosome dot. Require to shift manually the coordonates to add the chromosome labels not overlapping.
# annotations chr part 1
sub<-d[d$type=='sample',]
t<-data.frame(chr=1:22); t$y<-sub$rate[match(t$chr, sub$chr)]; t$x<-sub$cM[match(t$chr, sub$chr)]
t$y[t$chr%%2==0]<-t$y[t$chr%%2==0]+1
t$y[t$chr%%2!=0]<-t$y[t$chr%%2!=0]-1
t$y[t$chr==17]<-t$y[t$chr==17]+2
t$y[t$chr==11]<-t$y[t$chr==11]+2
t$y[t$chr==8]<-t$y[t$chr==8]-2
t$y[t$chr==9]<-t$y[t$chr==9]+2
t$y[t$chr==7]<-t$y[t$chr==7]+2
t$y[t$chr==6]<-t$y[t$chr==6]-2
t$y[t$chr==5]<-t$y[t$chr==5]+2
t$x[t$chr==18]<-t$x[t$chr==18]-6
t$y[t$chr==18]<-t$y[t$chr==18]-1
t$x[t$chr==20]<-t$x[t$chr==20]-6
t$y[t$chr==20]<-t$y[t$chr==20]-0.5
t$x[t$chr==22]<-t$x[t$chr==22]-8
t$y[t$chr==22]<-t$y[t$chr==22]-1
t$x[t$chr==16]<-t$x[t$chr==16]+1
t$x[t$chr==13]<-t$x[t$chr==13]+1
t$x[t$chr==8]<-t$x[t$chr==8]+1
t$x[t$chr==19]<-t$x[t$chr==19]+4
t$y[t$chr==19]<-t$y[t$chr==19]
t$x[t$chr==12]<-t$x[t$chr==12]-1
t$x[t$chr==11]<-t$x[t$chr==11]-1
t$x[t$chr==14]<-t$x[t$chr==14]-2
t$chr[t$chr==21]<-paste0('chr',t$chr[t$chr==21])
t$chr[t$chr==22]<-paste0('chr',t$chr[t$chr==22])
t$chr[t$chr==1]<-paste0('chr',t$chr[t$chr==1])
t$chr[t$chr==2]<-paste0('chr',t$chr[t$chr==2])
t$x[t$chr=='chr22']<-t$x[t$chr=='chr22']-3
t$y[t$chr=='chr22']<-t$y[t$chr=='chr22']+0.5
t1<-t

# annotations chr part 2
sub<-d[d$type=='call',]
t<-data.frame(chr=1:22); t$y<-sub$rate[match(t$chr, sub$chr)]; t$x<-sub$cM[match(t$chr, sub$chr)]
t$y[t$chr%%2==0]<-t$y[t$chr%%2==0]+1
t$y[t$chr%%2!=0]<-t$y[t$chr%%2!=0]-1
t$y[t$chr==7]<-t$y[t$chr==7]+2
t$y[t$chr==6]<-t$y[t$chr==6]-2
t$y[t$chr==10]<-t$y[t$chr==10]-2
t$y[t$chr==11]<-t$y[t$chr==11]+2
t$y[t$chr==13]<-t$y[t$chr==13]+2
t$y[t$chr==18]<-t$y[t$chr==18]-2
t$x[t$chr==8]<-t$x[t$chr==8]-2
t$x[t$chr==22]<-t$x[t$chr==22]-2
t$x[t$chr==14]<-t$x[t$chr==14]-4


t$chr[t$chr==21]<-paste0('chr',t$chr[t$chr==21])
t$chr[t$chr==22]<-paste0('chr',t$chr[t$chr==22])
t$chr[t$chr==1]<-paste0('chr',t$chr[t$chr==1])
t$chr[t$chr==2]<-paste0('chr',t$chr[t$chr==2])
t$x[t$chr=='chr22']<-t$x[t$chr=='chr22']-4
t$y[t$chr=='chr22']<-t$y[t$chr=='chr22']+0.25

t2<-t

gf<-ggplot(d, aes(x=cM, y=rate, col=type, group=type)) + geom_point(size=point_size-0.8) + geom_line(size=line_size-0.3) +
  scale_color_manual(values=c(mycols[4], mycols[5]), labels=c('% of call', '% of samples')) + theme_bw() + theme(legend.position = c(0.006,.994),legend.justification = c(0, 1),legend.direction="horizontal",
                                                                                                                 legend.box.background = element_rect(color="black", size=1),
                                                                                                                 panel.border = element_rect(color="black", size=1),
                                                                                                                 axis.ticks = element_line(size=.4), axis.text = element_text(size=st),
                                                                                                                 axis.title = element_text(size=st), legend.text = element_text(size=lt),
                                                                                                                 legend.title = element_text(size=lt), legend.key.height = unit(0,'cm'),
                                                                                                                 legend.key.width = unit(0,'cm'), legend.key.size = unit(0, 'cm'),
                                                                                                                 legend.margin = margin(5,5,5,5), legend.spacing.y = unit(0,'cm'),
                                                                                                                 title = element_text(size=stt),plot.title = element_text(face='bold',hjust = 0.5),
                                                                                                                 plot.margin = margin(1, 5.5, 5.5, 3.5),panel.grid = element_blank())+ 
  guides(col=guide_legend(nrow=1,byrow=TRUE, override.aes = list(size=.5))) +
  labs(y='Rate [%]', x='Length of chromosomes [cM]', col='') + lims(y=c(50,100)) +
  annotate("text", x=t1$x, y=t1$y, label=t1$chr, size=1.5) +
  annotate("text", x=t2$x, y=t2$y, label=t2$chr, size=1.5)



title <- ggdraw() + 
  draw_label('F. Rates per chromosome length',
             fontface = 'bold',hjust = 0.5  , size=title_size) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))



Figure_2F<-plot_grid(title, gf, ncol=1,rel_heights = c(0.03, 1))





assembly_DEF<-plot_grid( Figure_2D, Figure_2E, Figure_2F, ncol=3, rel_widths=c(1,1.1,1))



Figure2<-plot_grid(assembly_ABC2, assembly_DEF, ncol=1, rel_heights = c(1,1))

#ggsave(file="Final_figures/Figure2_Rplot2.png", assembly_fig2, width = a4_width, height = 150, units = "mm", dpi=600)






