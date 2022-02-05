tabletext <- cbind(c("Variable",hrtable$Variable),
                   c("HR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(hrtable$lower.95CI),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(hrtable$upper.95CI),3),nsmall = 3)),
                   c("Pvalue",formatC(as.numeric(hrtable$p), format = "e", digits = 2)))


tabletext[2,] <- c("Univariable Cox (George et.al.)",NA,NA,NA,NA)
tabletext[6,] <- c("Multivariable Cox (George et.al.)",NA,NA,NA,NA)


pdf("forestplot of risk table.pdf", width = 11, height = 6)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),#log2(HR)
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), 
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#a3a659", lines="black", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-1,0,1),
           lwd.xaxis=2,
           xlab=expression("log"[2]~"HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="black", lty=2),
                           "6" = gpar(lwd=1, col="black", lty=2),
                           "10" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.9,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())



gmtfile1<-'c2.cp.kegg.v7.1.symbols.gmt'
gmtfile2<-'c2.cp.reactome.v7.1.symbols.gmt'
gmtfile3<-'c5.all.v7.1.symbols.gmt'
d1 <- read.gmt(gmtfile1)
d2 <- read.gmt(gmtfile2)
d3 <- read.gmt(gmtfile3)
dd<-rbind(d1,d2,d3)
hallmark.list <- dd %>% split(.$ont) %>% lapply( "[[", 2)
save(hallmark.list,file = 'GSEA_List.Rdata')


gsym.fc<-read.csv("input.csv",
                  check.names = F)
gsym.fc<-gsym.fc[,c(2,3)]
ID <- bitr(gsym.fc$ENTREZID, fromType =c("ENTREZID"),
           toType =  "SYMBOL",
           OrgDb = org.Hs.eg.db) 

ID$ENTREZID<-as.numeric(ID$ENTREZID)
colnames(gsym.fc)
gsym.fc<-left_join(gsym.fc,ID,by=c("ENTREZID"))
gsym.fc<-gsym.fc[order(gsym.fc$logFC,decreasing = T),]

si.id <- gsym.fc$logFC; names(si.id) <- gsym.fc$SYMBOL

fgseaRes <- fgsea(pathways = hallmark.list, 
                  stats = si.id,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
pa<-read.csv('Final_Select_Pathways.csv',check.names = F)
ind<-read.csv('msigdb_v7.0/GSEA_ID.csv',check.names = F)
pa<-left_join(pa,ind,by=c('ID'))
pa<-pa[,-9]
pa<-na.omit(pa)
topPathways <- intersect(fgseaRes$pathway,pa$STANDARD_NAME)
topPathways<-topPathways[-4]
pdf("fgsea.pdf", width = 15, height = 8)
plotGseaTable(hallmark.list[topPathways], si.id, fgseaRes, 
              gseaParam = 0.3)
dev.off()

mycol <- rep(c("#CD919E","#CDB38B","#C1CDC1","#4F94CD","#9F79EE","#A9A9A9","#CDCD00","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)

#格式转换
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "Cohort")
dim(UCB_lodes)

ggplot(UCB_lodes,
       aes(x = x, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8) + 
  geom_stratum(alpha = .9,width = 1/10) + 
  geom_text(stat = "stratum", size = 3,color="black") + 
  
  
  scale_fill_manual(values = mycol) +
  
  xlab("") + 
  ylab("") +
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(fill = FALSE) 


tmp=data.frame(gene=data_plot[,2],group=data_plot[,1])
pv<-wilcox.test(tmp[which(tmp$group  == "Low"),"gene"],tmp[which(tmp$group == "High"),"gene"])$p.value
jco <- c("#4682B4","#548B54")

tmp$group<-as.character(tmp$group)
n1<-length(tmp$group[tmp$group=='Low'])
n2<-length(tmp$group[tmp$group=='High'])
tmp$group[tmp$group=='Low']<-paste0('Low (',n1,')')
tmp$group[tmp$group=='High']<-paste0('High (',n2,')')

pv.lab.2<-case_when(pv<0.05&pv>0.01~"*",pv<=0.01&pv>0.001~"**",
                    pv<=0.001&pv>0.0001~"***",
                    pv<=0.0001~"****",
                    pv>0.05~"ns")
lab<-pv.lab.2
label <- ggdraw() + draw_label(lab,size = 7,fontface = "bold")
plot_grid(label, 
          ggplot(data = tmp,aes(x =group , y = gene, fill =group))+
            scale_fill_manual(values = jco[2:1]) + 
            geom_violin(alpha=0.4, position = position_dodge(width = .75),
                        size=0.4, color="black") + 
            geom_boxplot(notch = TRUE, outlier.size = -1, 
                         color="black", lwd=0.4, alpha = 0.7)+ 
            geom_point(shape = 21, size=0.4, 
                       position = position_jitterdodge(), 
                       color="black", alpha=1)+ 
            theme_classic() +
            ylab('log2 (TMB+1)')+
            xlab('MMP9')  +
            theme(
              axis.ticks = element_line(size=0.1,color="black"),
              axis.ticks.length = unit(0.1,"cm"),
              legend.position = "none",
              axis.title = element_text(size = 7),
              axis.text.x = element_text(size = 7,angle = 15,colour = "black",hjust = 1,vjust = 1),
              axis.text.y = element_text(size = 7,color = "black")), 
          ncol=1, rel_heights=c(.2, 1))
