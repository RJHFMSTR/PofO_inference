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

#Figure 3A
d<-as.data.frame(data.table::fread('source_data_F3A.txt.gz', hea=T))


d$BP<-as.numeric(d$BP)
d$CHR<-as.numeric(d$CHR)
d$P_ADD<-as.numeric(d$P_ADD)
d$P_MAT<-as.numeric(d$P_MAT)
d$P_PAT<-as.numeric(d$P_PAT)
d$P_COMP<-as.numeric(d$P_COMP)
# modify coordinate to plot
s<-0;d$BPcum<-NA;nCHR<-length(unique(d$CHR))
nbp<-c();sig<-5e-08
for (i in unique(sort(d$CHR))){
  nbp[i]<-max(d$BP[d$CHR==i])
  d$BPcum[d$CHR==i]<-d$BP[d$CHR==i]+s
  s<-s+nbp[i]}

dd<-d
axis.set<- dd %>% group_by(CHR) %>% summarize(center= (max(BPcum)+min(BPcum))/2)
ylim<-abs(floor(log10(min(c(dd$P_ADD, dd$P_PAT, dd$P_MAT, dd$P_COMP))))) +2


g_add<-ggplot() +
  geom_point(data=dd,  aes(  x=BPcum, y=-log10(as.numeric(P_ADD)), color=as.factor(CHR)), size=0.6, alpha=0.75)+
  geom_point(data=dd[dd$SNP==SNP,],  aes(  x=BPcum, y=-log10(as.numeric(P_ADD))),fill='grey', col='black',size=1.6, alpha=0.75, shape=23)+
  scale_x_continuous(label=axis.set$CHR, breaks=axis.set$center) + scale_y_continuous(expand=c(0,0), limits = c(0, -log10(min(dd$P_ADD))+2)) +
  scale_size_continuous(range=c(0.5,3)) + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
  scale_color_manual(values=rep(c('grey26','black'), nCHR)) + theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks = element_line(size=.7), axis.text.x = element_blank(), axis.title.y = element_text(size=st),
        axis.text.y=element_text(size=st), axis.title.x = element_blank(),
        title = element_text(size=6),
        plot.title = element_text(face='bold',hjust = 0.5)) +
  labs(y = "-log10(Additive p-value)") + ggtitle('')




g_mat<-ggplot() +
  geom_point(data=dd,  aes(  x=BPcum, y=-log10(as.numeric(P_MAT)), color=as.factor(CHR)), size=0.6, alpha=0.75)+
  geom_point(data=dd[dd$SNP==SNP,],  aes(  x=BPcum, y=-log10(as.numeric(P_MAT))),fill='yellow', col='red3',size=1.6, alpha=0.75, shape=23)+
  scale_x_continuous(label=axis.set$CHR, breaks=axis.set$center) + scale_y_continuous(expand=c(0,0), limits = c(0, -log10(min(dd$P_MAT))+2)) +
  scale_size_continuous(range=c(0.5,3)) + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
  scale_color_manual(values=rep(c('rosybrown','red3'), nCHR)) + theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks = element_line(size=.7), axis.text.x = element_blank(), axis.title.y = element_text(size=st),
        axis.text.y=element_text(size=st), axis.title.x = element_blank(),
        title = element_text(size=stt),
        plot.title = element_text(face='bold',hjust = 0.5)) +
  labs(y = "-log10(Maternal p-value)")

g_pat<-ggplot() +
  geom_point(data=dd,  aes(  x=BPcum, y=-log10(as.numeric(P_PAT)), color=as.factor(CHR)), size=0.6, alpha=0.75)+
  geom_point(data=dd[dd$SNP==SNP,],  aes(  x=BPcum, y=-log10(as.numeric(P_PAT))),fill='yellow', col='royalblue4',size=1.6, alpha=0.75, shape=23)+
  scale_x_continuous(label=axis.set$CHR, breaks=axis.set$center) + scale_y_continuous(expand=c(0,0), limits = c(0, -log10(min(dd$P_PAT))+2)) +
  scale_size_continuous(range=c(0.5,3)) + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
  scale_color_manual(values=rep(c('royalblue4','navyblue'), nCHR)) + theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks = element_line(size=.7), axis.text.x = element_blank(), axis.title.y = element_text(size=st),
        axis.text.y=element_text(size=st), axis.title.x = element_blank(),
        title = element_text(size=stt),
        plot.title = element_text(face='bold',hjust = 0.5)) +
  labs(y = "-log10(Paternal p-value)")


#ylimmax=-log10(min(dd$P_COMP))+2
ylimmax=11
g_diff<-ggplot() +
  geom_point(data=dd,  aes(  x=BPcum, y=-log10(as.numeric(P_COMP)), color=as.factor(CHR)), size=0.6, alpha=0.75)+
  geom_point(data=dd[dd$SNP==SNP,],  aes(  x=BPcum, y=-log10(as.numeric(P_COMP))),fill='yellow', col='darkgreen',size=1.6, alpha=0.75, shape=23)+
  scale_x_continuous(label=axis.set$CHR, breaks=axis.set$center) +
  scale_y_continuous(expand=c(0,0), limits = c(0, ylimmax), breaks=c(0,5,10), labels=c(0,5,10)) + 
  scale_size_continuous(range=c(0.5,3)) + geom_hline(yintercept=-log10(sig), col='red', alpha=0.2) +
  scale_color_manual(values=rep(c('darkgreen','springgreen4'), nCHR)) + theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks = element_line(size=.7), axis.text.x = element_text(size=st), axis.title.y = element_text(size=st),
        axis.text.y=element_text(size=st), axis.title.x = element_text(size=st),
        title = element_text(size=stt),
        plot.title = element_text(face='bold',hjust = 0.5)) +
  labs(y = "-log10(Differential p-value)", x='Chromosomes')


title <- ggdraw() + 
  draw_label('Figure 4.',
             fontface = 'bold',x = 0,hjust = 0  , size=10) +  theme(plot.margin = margin(0, 0, 0, 7)  ) +
  theme(plot.background = element_rect(fill='white', colour = 'white'))


Figure_3A<-plot_grid(title,g_add, g_mat, g_pat, g_diff, ncol=1, rel_heights = c(0.08,1,1,1,1))

#ggsave(file=paste0("Final_figures/GWAS_scan_",pheno,".png"), assembly_A, width = a4_width, height = a4_height*0.6, units = "mm")





