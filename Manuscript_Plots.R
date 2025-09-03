dir <- '~/Partners HealthCare Dropbox/Ning Guan/plots-Alexa/6. CRISPR Screen'
library(readxl)
library(ggrepel)
library(ggplot2)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(EnhancedVolcano)
rm(list=ls())

drawscatter<-function(df1,df2,xlab='test1',ylab='test2',labels,figname='test.pdf'){
  srt_df2<-df2[match(df1$id,df2$id),]
  dat<-cbind.data.frame(gene=df1$id,v1=-log10(df1$neg.p.value),v2=-log10(srt_df2$neg.p.value))
  dat <- dat %>% mutate(label = ifelse(gene %in% labels,gene,NA))
  highlight_df <- dat %>% filter(gene %in% labels)
  p <- ggplot(dat, aes(v1, v2))+
    geom_point(color = "grey", size=1)+
    geom_point(data=highlight_df, aes(v1, v2), color='black', size=2)+
    xlim(0,max(dat$v1)+2)+
    ylim(0,max(dat$v2)+2)+
    xlab(xlab)+
    ylab(ylab)+
    geom_text_repel(
      aes(label = label),
      # family = "Poppins",
      size = 5,
      max.overlaps = 5,
      min.segment.length = 0, 
      seed = 42, 
      box.padding = 0.5,
      arrow = arrow(length = unit(0.010, "npc")),
      nudge_x = .5,
      nudge_y = .5,
    )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
  print(p)
}
newscatter<-function(df,xlab='test1',ylab='test2',label_genes,figname='test.pdf'){
  srt_df2<-df2[match(df$id,df2$id),]
  dat<-cbind.data.frame(gene=df$id,v1=srt_df2$neg.lfc,v2=df$`L2_ET3:1_lfc`)
  labels <- intersect(df$id,label_genes)
  labels <- labels[!(labels %in% c("Tmem165","Slc35a2"))]
  dat <- dat %>% mutate(label = ifelse(gene %in% labels,gene,NA))
  dat <- dat %>% mutate(labelt = ifelse(gene %in% c("Tmem165","Slc35a2"),gene,NA))
  highlight_df <- dat %>% dplyr::filter(gene %in% labels)
  highlight_df2 <- dat %>% dplyr::filter(gene %in% c("Tmem165","Slc35a2"))
  p <- ggplot(dat, aes(v1, v2))+
    geom_point(color = "grey", size=1)+
    geom_point(data=highlight_df, aes(v1, v2), color='blue', size=2)+
    geom_point(data=highlight_df2, aes(v1, v2), color='red', size=2)+
    xlim(min(dat$v1),0)+#xlim(0,4.6)+#xlim(0,max(dat$v1)+1)+#
    ylim(min(dat$v2),0)+#ylim(-0.5,max(dat$v2)+0.5)+#
    xlab(xlab)+
    ylab(ylab)+
    geom_text_repel(
      aes(label = label),
      # family = "Poppins",
      size = 5,
      max.overlaps = Inf,
      min.segment.length = 0, 
      seed = 42, 
      box.padding = 0.5,
      arrow = arrow(length = unit(0.010, "npc")),color="black")+
    geom_text_repel(
      aes(label = labelt),
      # family = "Poppins",
      size = 6,
      max.overlaps = Inf,
      min.segment.length = 0, 
      seed = 42, 
      box.padding = 0.5,
      nudge_x=-0.5,
      nudge_y=0.3,
      arrow = arrow(length = unit(0.010, "npc")),color="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size=15),axis.title = element_text(size=17))
          #panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  
  print(p)
}

mouse_gene_sets = msigdbr(species = "Mus musculus")
msigdbr_collections()
mouse_GO = mouse_gene_sets %>% dplyr::filter(gs_cat == "C5")
go_glycosylation_terms <- c("GOBP_GLYCOSYLATION", "GOBP_REGULATION_OF_PROTEIN_GLYCOSYLATION",
                            "GOBP_PROTEIN_N_LINKED_GLYCOSYLATION","GOMF_GLYCOSYLTRANSFERASE_ACTIVITY",
                            "GOMF_SIALYLTRANSFERASE_ACTIVITY")
glycosylation_table = mouse_GO %>% dplyr::filter(gs_name %in% go_glycosylation_terms)
glycosylation_genes = as.vector(unique(glycosylation_table$gene_symbol))

#PB3 ET1:4 and L2 ET3:1 **
df <- read_excel(file.path(dir,"minipool-in vitro/minipool_mergedlist_Tcellkilling.xlsx"))
newscatter(df, xlab=expression(paste("Log"[2]," Fold Change PB3")),y=expression(paste("-Log"[10]," (p value) PB3")),glycosylation_genes)

#PB3 ET1:4 and L2 ET3:1
top_pb3 = df[order(df$`PB3_ET1:4_p.value`),]$id[1:10]
top_L2 = df[order(df$`L2_ET3:1_p.value`),]$id[1:10]
newscatter(df, xlab=expression(paste("-Log"[10]," (p value) PB3")),y=expression(paste("-Log"[10]," (p value) 4T1-Mes")),union(top_pb3,top_L2))
newscatter(df, xlab=expression(paste("Log"[2]," Fold Change PB3")),ylab=expression(paste("Log"[2]," Fold Change 4T1-Mes")),union(top_pb3,top_L2))

#PB3 vs. PB2
top_pb3 = df[order(df$`PB3_ET1:4_p.value`),]$id[1:10]
top_pb2 = df[order(df$`PB2_ET1:7_p.value`),]$id[1:10]
newscatter(df, xlab=expression(paste("Log"[2]," Fold Change PB3")),y=expression(paste("Log"[2]," Fold Change PB2")),union(top_pb3,top_pb2))

# PB3 vs. 4T1 in vivo CTLA4 vs. NT
df2 <- read.csv(file.path(dir,"minipool-in vivo/rra_NT_CTLA4.4.gene_summary.txt"),sep="\t")
top_NT_CTLA4 <- df2[order(df2$neg.p.value),]$id[1:10]
newscatter(df, xlab=expression(paste("-Log"[10]," (p value) CTLA4 vs. NT")),y=expression(paste("-Log"[10]," (p value) 4T1-Mes")),union(top_L2,top_NT_CTLA4))
newscatter(df, xlab=expression(paste("Log"[2]," Fold Change CTLA4 vs. NT")),ylab=expression(paste("Log"[2]," Fold Change 4T1-Mes")),union(top_L2,top_NT_CTLA4))

df1<-read.csv('mini_pool-in vitro/PB3_ET_1:4_ET1_4_vsNoOVA_r1.gene_summary.txt',sep='\t')
df2<-read.csv('mini_pool-in vitro/L2_ET_3:1_ET3_1_vsNoGFP_r1.gene_summary.txt',sep='\t')
labels <- c("Tmem127", "Ap2s1", "Fcho2", "Susd6","Pigu", "Dpm2")
drawscatter(df1,df2,xlab='PB3_ET1:4',ylab='L2_ET3:1',labels, figname='figures/MiniPoolScreen/crossConditions/PB3_ET1:4vsL2ET3:1.pdf')

#PB3 ET1:4 and PB2_ET1:7 **
df1<-read.csv('mini_pool/PB3_ET_1:4_ET1_4_vsNoOVA_r1.gene_summary.txt',sep='\t')
df2<-read.csv('mini_pool/PB2_ET_1:7_ET1_7_vsNoOVA_r1.gene_summary.txt',sep='\t')
drawscatter(df1,df2,xlab='PB3_ET1:4',ylab='PB2_ET1:7',labels, figname='figures/MiniPoolScreen/crossConditions/PB3_ET1:4vsPB2_ET1:7.pdf')

# Make Bubble plot for KEGG
df1 <- read_excel(file.path(dir,"Primary screen/TNSC_Breast_OT-I.OVA_Coculture-T_r1+2+3.gene_summary.xlsx"),sheet="Xiaohan rep1")
df2 <- read_excel(file.path(dir,"Primary screen/TNSC_Breast_OT-I.OVA_Coculture-T_r1+2+3.gene_summary.xlsx"),sheet="Xiaohan rep2")
df3 <- read_excel(file.path(dir,"Primary screen/TNSC_Breast_OT-I.OVA_Coculture-T_r1+2+3.gene_summary.xlsx"),sheet="Xiaohan rep3")
srt_df2<-df2[match(df1$id,df2$id),]
srt_df3<-df3[match(df1$id,df3$id),]
dat<-cbind.data.frame(gene=df1$id,pval1=df1$`neg|p-value`,pval2=srt_df2$`neg|p-value`,pval3=srt_df3$`neg|p-value`,
                      lfc1 = df1$`neg|lfc`, lfc2=srt_df2$`neg|lfc`, lfc3=srt_df3$`neg|lfc`, 
                      posval1 = df1$`pos|p-value`,posval2 = df2$`pos|p-value`,posval3=df3$`pos|p-value`)
dat$med.lfc = rowMeans(dat[c('lfc1', 'lfc2','lfc3')])
neg_genes <- dat[dat$med.lfc < 0 & ((dat$pval1 < 0.05 & dat$pval2 <0.05)|
                                      (dat$pval1 < 0.05 & dat$pval3 < 0.05)|
                                      (dat$pval2 < 0.05 & dat$pval3 < 0.05)),]
pos_genes <- dat[dat$med.lfc > 0 & ((dat$posval1 < 0.05 & dat$posval2 <0.05)|
                                      (dat$posval1 < 0.05 & dat$posval3 < 0.05)|
                                      (dat$posval2 < 0.05 & dat$posval3 < 0.05)),]
# GO analysis using ClusterProfiler
matchMMUsymbol<-function(mouse_symbol){
  dir = "~/Partners HealthCare Dropbox/Ning Guan/Documents/Code"
  map=read.table(file.path(dir,'codes_Davis/src/Mouse_Human_GeneIDs_sub.txt'),sep=' ',header=T)
  map_sor<-map[match(mouse_symbol,map$mmu_symbol),]
  human<-map_sor$hsa_symbol
  return(human)
}

ngenes <- neg_genes$gene
human_ngenes <- matchMMUsymbol(ngenes)
human_symbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = human_ngenes,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "SYMBOL")
human_symbol <- human_symbol[!is.na(human_symbol$ENTREZID),]
name_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
       keys = ngenes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
name_symbol <- name_symbol[!is.na(name_symbol$ENTREZID),]
write.csv(ngenes,file.path(dir,"ngenes.csv"),quote=F,col.names = F,row.names = F)
x <- enrichKEGG(name_symbol$ENTREZID, organism="mmu")
human_kegg <- enrichKEGG(human_symbol$ENTREZID, organism="hsa")
# Make BubblePlot
x@result$Description <- gsub('.{29}$', "",x@result$Description)
x@result$Log10FDR <- -log10(x@result$p.adjust)
dotplot(x, x="Log10FDR",showCategory=10,font.size=11) + ggtitle("Genome Screen Negative Hits GO KEGG")
# Get geneID and convert to names
pgenes <- pos_genes$gene
pos_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keys = pgenes,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
pos_kegg <- enrichKEGG(pos_symbol$ENTREZID, organism="mmu")
pos_kegg@result$Description <- gsub('.{29}$', "",pos_kegg@result$Description)
pos_kegg@result$Log10FDR <- -log10(pos_kegg@result$p.adjust)

# Make GO plot for 4T1 Tmem165-KO vs. Ctrl
tmem165_ctrl <- read_excel("~/Partners HealthCare Dropbox/Ning Guan/plots-Alexa/
                           7. 4T1-Tmem165 vs Ctrl bulk RNA-seq/Tmem165_vs_Int2.deseq2_significant.xlsx")
tmem165_ctrl <- Tmem165_vs_Int2_deseq2_significant
colnames(tmem165_ctrl)[1] <- "Gene"
# Upregulated vs. Downregulated genes
ngenes <- tmem165_ctrl[tmem165_ctrl$log2FoldChange<0 & tmem165_ctrl$padj<0.05,]$Gene
pgenes <- tmem165_ctrl[tmem165_ctrl$log2FoldChange>0 & tmem165_ctrl$padj<0.05,]$Gene
# Bubbleplot
neg_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keys = ngenes,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
neg_symbol <- neg_symbol[!is.na(neg_symbol$ENTREZID),]
# Hallmark Gene Set
hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
neg_enriched <- enricher(gene = neg_symbol$ENTREZID,
                     TERM2GENE = hallmark,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
neg_enriched@result$Log10FDR <- -log10(neg_enriched@result$p.adjust)
dotplot(neg_enriched, x="Log10FDR",color="p.adjust",showCategory=10,font.size=11) + 
  ggtitle("Ctrl Upregulated GO Hallmark")
# Positive
pos_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                    keys = pgenes,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
pos_symbol <- pos_symbol[!is.na(pos_symbol$ENTREZID),]
pos_enriched <- enricher(gene = pos_symbol$ENTREZID,
                         TERM2GENE = hallmark,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)
pos_enriched@result$Log10FDR <- -log10(pos_enriched@result$p.adjust)
dotplot(pos_enriched, x="Log10FDR",color="p.adjust",showCategory=10,font.size=11) + 
  ggtitle("Tmem165 Upregulated GO Hallmark")


# Primary screen
merge_df <- read_excel(file.path(dir,"Primary screen/TNSC_Breast_OT-I.OVA_Coculture-T_r1+2+3.gene_summary.xlsx"),
                  sheet="Xiaohan_TNBC_Breast_OT-I.OVA_me")
merge_df <- transform(merge_df, pvalue = pmin(pos.p.value, neg.p.value))
merge_df$id[5] <- "Mosmo"
merge_df <- transform(merge_df, fdr = pmin(pos.fdr, neg.fdr))
# Get genes in each pathway
# Merge x and pos_kegg
merged_kegg <- rbind(head(as.data.frame(x@result)),head(as.data.frame(pos_kegg@result)))
# Convert geneID to vector of gene names
merged_kegg$genes <- strsplit(merged_kegg$geneID,"\\/")
df <- data_frame(matrix(nrow = 0, ncol = 3))
colnames(df) <- c("ENTREZID","SYMBOL","Pathway")
for (i in 1:12){
  geneid <- as.vector(unlist(merged_kegg$genes[i]))
  gene_names <- AnnotationDbi::select(org.Mm.eg.db, 
                                 keys = geneid,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "ENTREZID")
  gene_names$Pathway <- merged_kegg$Description[i]
  df <- rbind(df,gene_names)
}
ts <- AnnotationDbi::select(org.Mm.eg.db, 
                                          keys = c("21982","22232"),
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "ENTREZID")
ts$Pathway <- "N-Glycan biosynthesis"
df <- rbind(df,ts)# Add Tmem165 and Slc35a2
srt_df <-df[match(merge_df$id,df$SYMBOL),]
#srt_df$Pathway <- factor(srt_df$Pathway,levels=unique(df$Pathway))
dat<-cbind.data.frame(gene=merge_df$id,v1=merge_df$neg.lfc,v2=-log10(merge_df$pvalue),
                      pathway=srt_df$Pathway,pvalue = merge_df$pvalue)
dat <- dat %>% mutate(pathway = ifelse(pvalue < 0.01|gene=="Tmem165",pathway,NA))
dat$pathway <- factor(dat$pathway, levels = c("N-Glycan biosynthesis", 
                  "Protein processing in endoplasmic reticulum", "Autophagy - animal",
                  "Terpenoid backbone biosynthesis", "Endocytosis",
                 "Human T-cell leukemia virus 1 infection",
                  "Basal transcription factors",NA))
dat <- dat %>% mutate(label = ifelse((gene %in% df$SYMBOL & pvalue < 0.01)|(gene %in% c("Tmem165","Slc35a2")),gene,NA))
highlight_df <- dat %>% dplyr::filter((gene %in% df$SYMBOL & pvalue < 0.01)|(gene %in% c("Tmem165","Slc35a2")))
p <- ggplot(dat, aes(v1, v2))+
  geom_point(color = "grey", size=1)+
  geom_point(data=highlight_df, aes(v1, v2,col=pathway), size=1)+
  xlim(min(dat$v1),max(dat$v1))+#xlim(0,4.6)+#xlim(0,max(dat$v1)+1)+#
  ylim(0,max(dat$v2))+#ylim(-0.5,max(dat$v2)+0.5)+#
  labs(title = "PB3 Genome-scale Screen",x = expression(paste("Log"[2]," Fold Change")),
       y = expression(paste("-Log"[10]," (p value)")))+
  geom_vline(xintercept=0,color="grey50")+
  geom_text_repel(
    aes(label = label,col=pathway),
    # family = "Poppins",
    size = 5,
    max.overlaps = Inf,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    arrow = arrow(length = unit(0.010, "npc")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=15),axis.title = element_text(size=17),
        plot.title = element_text(hjust=0.5,size=17),legend.text = element_text(size=13),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position = "none")

EnhancedVolcano(dat,lab=dat$id,x="pos.lfc", y="pvalue",labSize=4,
                title="PB3 Genome-Scale Screen",pointSize=1,drawConnectors = T)

EnhancedVolcano(crispr1_neg, rownames(crispr1_neg), 'neg.lfc', 'neg.p.value', pCutoff = 6.5e-05,
                FCcutoff=0.5, labSize=4,
                pointSize = c(ifelse(crispr1_neg$neg.p.value<6.5e-05, 4, 1)),
                xlim = c(-3,0), ylim = c(0,7),
                col = c("grey30", "forestgreen", "royalblue", "red2"), colAlpha = 1/2,
                colGradientBreaks = c(pCutoff, 1), colGradientLimits = c(0, 1),
                axisLabSize = 15, titleLabSize = 20, subtitleLabSize = 8, captionLabSize = 8,
                legendPosition='none',
                title = 'Genes with depleted gRNA', subtitle=NULL,
                gridlines.major = FALSE, gridlines.minor = FALSE,
                border = "full",
                borderWidth = 0.8,
                borderColour = "black",
                drawConnectors = TRUE)

# Make Scatter Plot of Spearman correlation between Tmem165, CD8+ T, Tregs, MDSCs, and DCs
# Use Cybersort first, then TIDE otherwise. 
dir1 <- "~/Partners HealthCare Dropbox/Ning Guan/Lab Work- Alexa Mitchell/Spring_2025/TCGA BRCA Analysis/TCGA_Gene_Immune_Correlation"
tmem165_cd8 <- read.csv(file.path(dir1, "Tmem165_CD8_gene_table.csv"))
tmem165_treg_mdsc <- read.csv(file.path(dir1, "Tmem165_Tregs_MDSC_gene_table.csv"))
slc35a2_cd8 <- read.csv(file.path(dir1, "Slc35a2_CD8_gene_table.csv"))
Slc35a2_treg_mdsc <- read.csv(file.path(dir1, "Slc35a2_Tregs_MDSC_gene_table.csv"))
tmem165_cd4 <-  read.csv(file.path(dir1, "Tmem165_CD4_table.csv"))
tmem165_nk <-  read.csv(file.path(dir1, "Tmem165_NK_table.csv"))
tmem165_dc <-  read.csv(file.path(dir1, "Tmem165_DC_table.csv"))
slc35a2_cd4 <- read.csv(file.path(dir1, "Slc35a2_CD4_table.csv"))
slc35a2_nk <- read.csv(file.path(dir1, "Slc35a2_NK_table.csv"))
slc35a2_dc <- read.csv(file.path(dir1, "Slc35a2_DC_table.csv"))

# Tmem165 correlation with immune infiltration 
ndf <- tmem165_cd8
tcga_scatter <- function(ndf, xlab, ylab,title){
ndf$cancer_type <- trimws(gsub("\\(.*","",ndf$cancer))
# Select CyberSort, remove remaining rows
dat <- ndf[grepl("activated_CIBERSORT", ndf$infiltrates, fixed = TRUE),]
dat <- dat[!grepl("-ABS", dat$infiltrates, fixed = TRUE),]
dat$v1 = (-log10(dat$p))
dat$v2 = (-log10(dat$adj.p))
labels <- dat[dat$adj.p < 0.05,]$cancer_type
labels <- union(labels, "BRCA")
dat <- dat %>% mutate(label = ifelse(cancer_type %in% labels,cancer_type,NA))
highlight_df <- dat %>% dplyr::filter(cancer_type %in% labels)
p <- ggplot(dat, aes(rho, v1))+
  geom_point(color = "grey", size=1)+
  geom_point(data=highlight_df, aes(rho, v1), color='red', size=1)+
  xlim(min(dat$rho),max(dat$rho))+#xlim(0,4.6)+#xlim(0,max(dat$v1)+1)+#
  ylim(0,max(dat$v1))+#ylim(-0.5,max(dat$v2)+0.5)+#
  labs(title = title,x = xlab, y = ylab)+
  geom_text_repel(
    aes(label = label),
    # family = "Poppins",
    size = 5,
    max.overlaps = 30,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5)+
    #arrow = arrow(length = unit(0.010, "npc")),color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=15),axis.title = element_text(size=17),
        plot.title = element_text(hjust=0.5,size=17),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
print(p)
}
tcga_scatter(tmem165_cd8,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Tmem165 vs. CD8 T Cell") #CIBERSORT
tcga_scatter(tmem165_treg_mdsc,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Tmem165 vs. MDSC") #TIDE
tcga_scatter(slc35a2_cd8,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. CD8 T Cell")
tcga_scatter(Slc35a2_treg_mdsc,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. MDSC")
tcga_scatter(Slc35a2_treg_mdsc,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. Treg") #CIBERSORT
tcga_scatter(tmem165_cd4,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Tmem165 vs. CD4 T Cell") # QUANTISEQ
tcga_scatter(slc35a2_cd4,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. CD4 T Cell")
tcga_scatter(tmem165_nk,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Tmem165 vs. Activated NK Cell") #CIBERSORT
tcga_scatter(slc35a2_nk,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. Activated NK Cell")
tcga_scatter(tmem165_dc,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Tmem165 vs. Activated DC") #CIBERSORT
tcga_scatter(slc35a2_dc,"Spearman Correlation",
             expression(paste("-Log"[10]," (p value)")), "Slc35a2 vs. Activated DC")
