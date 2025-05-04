if (!require("BiocManager", quietly = TRUE))
{
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
  BiocManager::install("edgeR")
  BiocManager::install("EnhancedVolcano")
  BiocManager::install("M3C")
  BiocManager::install("ComplexHeatmap")
  BiocManager::install("fgsea")
  BiocManager::install("biocLite")
  n#library(org.Hs.eg.db)
}
library(svglite)
library(umap)
library(Rtsne)

library(edgeR)
library(survminer)
library(biomaRt)
library(circlize)
library(ComplexHeatmap)
library(DESeq2)
library(EnhancedVolcano)
library(TCGAbiolinks)
library(tidyverse)
library(survival)
library(dplyr)
library(tidyr)
#library(DTUrtle)
#library(tibble)
# #remotes::install_github("TobiTekath/DTUrtle"), install git from its website, install dependencies from vignette
# library(msigdbr)
# library(ClassDiscovery)
# library(Seurat)
# #remotes::install_github("mojaveazure/seurat-object", "seurat5")
# library(BisqueRNA)
# library(VennDiagram)
 library(magick)
library(cowplot)
#library("ggpubr")
# library(grid)
 library(gridExtra)
# library(M3C)
# library (pheatmap)
# library(RColorBrewer)
# library(LW1949)
 library(stringr)
 library(readr)
library(rstatix)
# library(rgl)

tricolor <- c("#4575B4", "#FFFFBF", "#D73027")

sigres <- function(res, padj_l = 0.01, l2fc_l = 1, indx=0) { #return significant result by p_limit and log_limit
  ##Ex ## sigres(res, 0.01, 2, 1) 
  if (indx == 0) { 
    subset(res, padj<padj_l & abs(log2FoldChange)>l2fc_l)
  } else {
    which (res$padj < padj_l)[which(abs(res[which (res$padj < padj_l),][,"log2FoldChange"]) > l2fc_l)]
  }
} #return list of significant results by threshold

EVplot <- function (myres){
  myVP<-EnhancedVolcano(myres,
                        lab = rownames(myres),
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.01,
                        FCcutoff = 2,
                        labSize = 5.0,
                        pointSize = 3.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.0)
  dev.new()
  print(myVP)
} #volcano plot for results(dds) of deseq

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = TRUE)
  breaks[!duplicated(breaks)]
}

quantile_breaks2 <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(min(xs), max(xs), length.out = n), na.rm = TRUE)
  breaks[!duplicated(breaks)]
}

scale_rows <- function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

# a long path with comma between vlaues
rmats_path_generator <- function (samples, pre, post){
  paste(paste0(pre,samples,post), collapse=',')
}

pcaMethods2 <- function (object, # deseq object
                        objres,  #deseq result object
                        intgroup="State", # by which column to divide groups
                        ntop=500, #number of genes for pca
                        returnData=FALSE, 
                        method = "var", #by which result perfrom pca
                        list = c(1), #list of possible genes
                        pc1 = 1, pc2 = 2, #pcs to use
                        pCutoff = 0.01, #for deseq res genes
                        FCcutoff = 1, #for deseq res genes
                        pSize = 1, #point size
                        labStyle = 1, #label type
                        colorStyle = 1, #0 for normal size 1 for size according to inclusion level & colored as well
                        toScale = FALSE, #scale pca
                        lSize = 5, #size of label
                        tSize = 15,#size of label text
                        gStyle = 1, # 0 none 1 scaling colors
                        incLevel_size = 1, #for feature plot need vector of incLevels
                        sStyle = 0, # 0 none 1 with chagning shapes
                        s2Style = FALSE, #2 true for 1 shape feature plot, choose shape number (21 - circle with fill)
                        lColor = 1, #color of label 1 is rainbow else is black
                        marking1 = FALSE, #mark selected points1
                        marking2 = FALSE,#mark selected points2
                        myTitle = "Graph",
                        color_wheel = tricolor,#c("#4575B4","#FFFFBF", "#D73027"),
                        groupMeta = NULL){
  if (length(incLevel_size) > 1)
    incLevel_size <- (incLevel_size - min(incLevel_size)) / (max(incLevel_size) - min(incLevel_size))
  
  # calculate the variance for each gene
  switch (method,
          "all" = filter <- 1:length(counts(object, normalized = TRUE)[,1]),
          "var" = filter <- rowVars(counts(object, normalized = TRUE)), # object is dds
          "signifi" = filter <- sigres (objres,pCutoff, FCcutoff, 1 ),
          "list" = filter <- list) # object is res
  
  # select the ntop genes
  select <- order(filter, decreasing=TRUE)[seq_len(min(ntop, length(filter)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  df_pca <- t(counts(object, normalized = TRUE)[select,])
  if (toScale == FALSE) pca <- prcomp(df_pca)
  if (toScale == TRUE) pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  if (!is.null(groupMeta))
    group <-factor(patient_state(groupMeta), 
                   levels = c('Normal', 
                              'Cancer', 
                              'Primary Tumor', 
                              'Thrombus Tumor'))
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2], group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ")
  ys = paste0("PC", as.character(pc2),": ")
  
  #create color for lables
  if (lColor == 1)
    lColor <- rainbow(length(levels(d$group)))[as.integer(d$group)]
  else
    lColor <- "black"
  
  cbreaks <- quantile_breaks(incLevel_size,101)
  switch (labStyle,
          lStyle <- geom_blank(),
          lStyle <- geom_text(label=colnames(counts(object, normalized = TRUE)),hjust=0, vjust=0, size = lSize,color = lColor),
          lStyle <- geom_label(size = lSize, label=colnames(counts(object, normalized = TRUE)), nudge_x = 0.25, nudge_y = 0.25, check_overlap = T),
          lStyle <- geom_label_repel(size = lSize,max.overlaps = Inf,aes(label = colnames(counts(object, normalized = TRUE))), box.padding   = 0.1, point.padding = 0.1,segment.color = 'grey50'))
  switch (colorStyle,
          cStyle <- geom_point(size=pSize * incLevel_size),
          cStyle <- geom_point(aes(color = incLevel_size), size=pSize * log(100*(incLevel_size+0.0101))),
          cStyle <- geom_point(aes(fill = incLevel_size), colour = 'black', size=pSize * log(100*(incLevel_size+0.0101)))) 
  
  switch (gStyle,
          grStyle <- geom_blank(),
          grStyle <- scale_colour_gradientn(colors = colorRampPalette(color_wheel)(length(cbreaks)-1),
                                            name = "Inclusion level",
                                            values = cbreaks),
          grStyle <- scale_fill_gradientn(colors = colorRampPalette(color_wheel)(length(cbreaks)-1),
                                            name = "Inclusion level",
                                            values = cbreaks))
  if(class(marking1) == "logical")
  {
    marking1 <- geom_blank()
  }
  else
  {
    marking1 <- geom_point(data = d[marking1,],
                          aes(x = PC1, y= PC2),
                          pch=21, 
                          fill=NA, 
                          size=pSize+3, colour="darkorange", stroke=0.5)
  }
  if(class(marking2) == "logical")
  {
    marking2 <- geom_blank()
  }
  else
  {
    marking2 <- geom_point(data = d[marking2,],
                          aes(x = PC1, y= PC2),
                          pch=21, 
                          fill=NA, 
                          size=pSize+3, colour="deepskyblue", stroke=0.5)
  }
  if ((length(unique(d$group)) < 9) & (sStyle == 1))
  {
    shape_values <- c(15,16,17,11,13,8,9,12)[1:length(unique(d$group))]
    tguide <- "legend"
    if(s2Style!= FALSE) {
      shape_values <- rep(s2Style,length(unique(d$group)))
      tguide <- "none"}
    sStyle <- scale_shape_manual(values = shape_values, guide = tguide)
    ggplot(data=d, aes(x=d$PC1, y=d$PC2, shape = d$group)) + 
      lStyle +
      cStyle + 
      grStyle + 
      sStyle + 
      marking1 +
      marking2 +
      xlab(paste0(xs,round(percentVar[pc1] * 100),"% variance")) +
      ylab(paste0(ys,round(percentVar[pc2] * 100),"% variance")) +
      labs(shape = paste(intgroup, collapse = ":")) + 
      coord_fixed() + ggtitle(myTitle) + theme(plot.title = element_text(size=tSize),
                                               legend.title=element_text(size=tSize))
  }
  else
  {
    legend_title <- ifelse(!is.null(groupMeta), 'State', paste(intgroup, collapse = ":"))
    sStyle <-geom_blank()
  ggplot(data=d, aes(x=d$PC1, y=d$PC2, color = d$group)) + 
    labs(type = paste(intgroup, collapse = ":")) + #glossary for color group
    lStyle +
    cStyle + 
    grStyle + 
    sStyle + 
    marking1 +
    marking2 +
    xlab(paste0(xs,round(percentVar[pc1] * 100),"% variance")) +
    ylab(paste0(ys,round(percentVar[pc2] * 100),"% variance")) +
    labs(color = legend_title) + 
    coord_fixed() + ggtitle(myTitle) + theme(plot.title = element_text(size=tSize),
                                             legend.title=element_text(size=tSize))
  }
  
}
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#data of the form 1-srr,2-name
igv_session_gen <- function(n,c,w,p,t,loc="All",gene="NONE", dir_path = "NONE", path_end = "")
{
  session_name <- paste0(dir_path,"/igv_",n[2],'_',c[2],'_',w[2],'_',p[2],'_',t[2],'_',gene,path_end,".xml")#gsub(':','-',format(Sys.time(), "%X"))
  session_code <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
<Session genome=\"hg38\" locus=\"",loc,"\" version=\"8\">
    <Resources>
        <Resource index=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",n[1],"/Aligned.sortedByCoord.out.bam.bai\" path=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",n[1],"/Aligned.sortedByCoord.out.bam\" type=\"bam\"/>
        <Resource index=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",c[1],"/Aligned.sortedByCoord.out.bam.bai\" path=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",c[1],"/Aligned.sortedByCoord.out.bam\" type=\"bam\"/>
        <Resource index=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",w[1],"/Aligned.sortedByCoord.out.bam.bai\" path=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",w[1],"/Aligned.sortedByCoord.out.bam\" type=\"bam\"/>
        <Resource index=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",p[1],"/Aligned.sortedByCoord.out.bam.bai\" path=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",p[1],"/Aligned.sortedByCoord.out.bam\" type=\"bam\"/>
        <Resource index=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",t[1],"/Aligned.sortedByCoord.out.bam.bai\" path=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",t[1],"/Aligned.sortedByCoord.out.bam\" type=\"bam\"/>
    </Resources>
    <Panel height=\"161\" name=\"Panel1684230233159\" width=\"1322\">
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Coverage\" autoScale=\"true\" clazz=\"org.broad.igv.sam.CoverageTrack\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",n[1],"/Aligned.sortedByCoord.out.bam_coverage\" name=\"Aligned.sortedByCoord.out.bam Coverage\" snpThreshold=\"0.2\" visible=\"true\">
            <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"60.0\" minimum=\"0.0\" type=\"LINEAR\"/>
        </Track>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Junctions\" autoScale=\"false\" clazz=\"org.broad.igv.sam.SpliceJunctionTrack\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",n[1],"/Aligned.sortedByCoord.out.bam_junctions\" maxdepth=\"50\" name=\"Aligned.sortedByCoord.out.bam Junctions\" visible=\"true\"/>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam\" clazz=\"org.broad.igv.sam.AlignmentTrack\" color=\"185,185,185\" displayMode=\"EXPANDED\" experimentType=\"RNA\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",n[1],"/Aligned.sortedByCoord.out.bam\" name=\"",n[2],"\" visible=\"true\">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height=\"161\" name=\"Panel1684230262375\" width=\"1322\">
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Coverage\" autoScale=\"true\" clazz=\"org.broad.igv.sam.CoverageTrack\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",c[1],"/Aligned.sortedByCoord.out.bam_coverage\" name=\"Aligned.sortedByCoord.out.bam Coverage\" snpThreshold=\"0.2\" visible=\"true\">
            <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"60.0\" minimum=\"0.0\" type=\"LINEAR\"/>
        </Track>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Junctions\" autoScale=\"false\" clazz=\"org.broad.igv.sam.SpliceJunctionTrack\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",c[1],"/Aligned.sortedByCoord.out.bam_junctions\" maxdepth=\"50\" name=\"Aligned.sortedByCoord.out.bam Junctions\" visible=\"true\"/>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam\" clazz=\"org.broad.igv.sam.AlignmentTrack\" color=\"185,185,185\" displayMode=\"EXPANDED\" experimentType=\"RNA\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",c[1],"/Aligned.sortedByCoord.out.bam\" name=\"",c[2],"\" visible=\"true\">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height=\"161\" name=\"Panel1684230295584\" width=\"1322\">
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Coverage\" autoScale=\"true\" clazz=\"org.broad.igv.sam.CoverageTrack\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",w[1],"/Aligned.sortedByCoord.out.bam_coverage\" name=\"Aligned.sortedByCoord.out.bam Coverage\" snpThreshold=\"0.2\" visible=\"true\">
            <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"60.0\" minimum=\"0.0\" type=\"LINEAR\"/>
        </Track>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Junctions\" autoScale=\"false\" clazz=\"org.broad.igv.sam.SpliceJunctionTrack\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",w[1],"/Aligned.sortedByCoord.out.bam_junctions\" maxdepth=\"50\" name=\"Aligned.sortedByCoord.out.bam Junctions\" visible=\"true\"/>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam\" clazz=\"org.broad.igv.sam.AlignmentTrack\" color=\"185,185,185\" displayMode=\"EXPANDED\" experimentType=\"RNA\" fontSize=\"10\" id=\"/home/be/zeevc/papers/ccRCC_Wang_2020/alignCa/",w[1],"/Aligned.sortedByCoord.out.bam\" name=\"",w[2],"\" visible=\"true\">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height=\"161\" name=\"Panel1684230346121\" width=\"1322\">
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Coverage\" autoScale=\"true\" clazz=\"org.broad.igv.sam.CoverageTrack\" fontSize=\"10\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",t[1],"/Aligned.sortedByCoord.out.bam_coverage\" name=\"Aligned.sortedByCoord.out.bam Coverage\" snpThreshold=\"0.2\" visible=\"true\">
            <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"60.0\" minimum=\"0.0\" type=\"LINEAR\"/>
        </Track>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Junctions\" autoScale=\"false\" clazz=\"org.broad.igv.sam.SpliceJunctionTrack\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",t[1],"/Aligned.sortedByCoord.out.bam_junctions\" maxdepth=\"50\" name=\"Aligned.sortedByCoord.out.bam Junctions\" visible=\"true\"/>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam\" clazz=\"org.broad.igv.sam.AlignmentTrack\" color=\"185,185,185\" displayMode=\"EXPANDED\" experimentType=\"RNA\" fontSize=\"10\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",t[1],"/Aligned.sortedByCoord.out.bam\" name=\"",p[2],"\" visible=\"true\">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height=\"161\" name=\"Panel1684230373268\" width=\"1322\">
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Coverage\" autoScale=\"true\" clazz=\"org.broad.igv.sam.CoverageTrack\" fontSize=\"10\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",p[1],"/Aligned.sortedByCoord.out.bam_coverage\" name=\"Aligned.sortedByCoord.out.bam Coverage\" snpThreshold=\"0.2\" visible=\"true\">
            <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"60.0\" minimum=\"0.0\" type=\"LINEAR\"/>
        </Track>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam Junctions\" autoScale=\"false\" clazz=\"org.broad.igv.sam.SpliceJunctionTrack\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",p[1],"/Aligned.sortedByCoord.out.bam_junctions\" maxdepth=\"50\" name=\"Aligned.sortedByCoord.out.bam Junctions\" visible=\"true\"/>
        <Track attributeKey=\"Aligned.sortedByCoord.out.bam\" clazz=\"org.broad.igv.sam.AlignmentTrack\" color=\"185,185,185\" displayMode=\"EXPANDED\" experimentType=\"RNA\" fontSize=\"10\" id=\"/scratch/zeevc/papers/ccRCC_Wang_2020/aligntk/",p[1],"/Aligned.sortedByCoord.out.bam\" name=\"",t[2],"\" visible=\"true\">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height=\"40\" name=\"FeaturePanel\" width=\"1322\">
        <Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" shouldShowTranslation=\"false\" visible=\"true\"/>
        <Track attributeKey=\"Refseq Genes\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;836.0;255,255,255;0,0,178\" fontSize=\"10\" groupByStrand=\"false\" id=\"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz\" name=\"Refseq Genes\" visible=\"true\"/>
    </Panel>
    <PanelLayout dividerFractions=\"0.007444168734491315,0.18982630272952852,0.37965260545905705,0.5694789081885856,0.7593052109181141,0.9491315136476427\"/>
    <HiddenAttributes>
        <Attribute name=\"DATA FILE\"/>
        <Attribute name=\"DATA TYPE\"/>
        <Attribute name=\"NAME\"/>
    </HiddenAttributes>
</Session>")
  
  if(!dir.exists(dir_path))dir.create(dir_path)
  writeLines(session_code, session_name) 
}

patient_state <- function (meta)
{
  (ifelse(meta$Normal == "Yes",
          "Normal",(ifelse(meta$Thrombosis == "No", 
                           "Cancer", ifelse(meta$State == "Cancer",
                                            "Primary Tumor","Thrombus Tumor")))))
}
dist2 <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

color_hm_ann <- function(df,hm,rl1 = NULL, rl2 = NULL, cl1 = NULL, cl2 = NULL, h = 70, w = 20)
{
  df$colors <- ifelse(rownames(df) %in% rl1,"orangered",
                         ifelse(rownames(df) %in% rl2,"royalblue","gray"))
  cols=df[order(match(rownames(df), hm$gtable$grobs[[5]]$label)), ]$colors
  hm$gtable$grobs[[6]]$gp=gpar(col=cols)
  save_pheatmap_pdf(hm, paste0('plots/hm','.pdf'), width = w,height = h)
}

ggSashimi_loc_gen <- function(events_path, out_name = "Loc_list.txt", file_path = T)
{
  if (file_path){rmats_events <- read.csv(events_path, row.names = 1)
  }else {rmats_events <- events_path}
  event_names <- rownames(rmats_events)
  ggS_event_name <- event_names#sapply(event_names,function (name) substr(name,1, unlist(gregexpr('_',name, fixed = T))[5]-1))
  events_loc <- paste0(rmats_events[,"chr"],':',
                rmats_events[,"upstreamES"],'-',
                rmats_events[,"downstreamEE"],'\n',ggS_event_name)
  events_loc_string <- paste0(events_loc,collapse = '\n')
  cat(events_loc_string, file = out_name)
}

ggSashimi_loc_rmats <- function(events_path, out_name = "Loc_list.txt")
{
  rmats_events <- read.csv(events_path, header = T)
  ggS_event_name <- apply(rmats_events,1, function (name) paste0(name['geneSymbol'],'_',
                                                                 name['chr'],'_',
                                                                 name['strand'],'_',
                                                                 name['exonStart_0base'],'_',
                                                                 name['exonEnd']))
  ggS_event_name <- gsub(' ', '', ggS_event_name, fixed = T)
  events_loc <- paste0(rmats_events[,"chr"],':',
                       rmats_events[,"upstreamES"],'-',
                       rmats_events[,"downstreamEE"],'\n',ggS_event_name)
  events_loc_string <- paste0(events_loc,collapse = '\n')
  cat(events_loc_string, file = out_name)
}

gene_lists_gen <- function(ctstk, wrs_res)
{
  #gene lists for sub-analysis
  mesenchymal_promoter_genes <- read.csv("gene_lists/kidney_emt.csv")
  mesenchymal_promoter_genes <- mesenchymal_promoter_genes[mesenchymal_promoter_genes$HUGO.Gene.Symbol %in% rownames(ctstk),]$HUGO.Gene.Symbol
  motif_genes <- read.csv("gene_lists/enriched motifs.Ver-4.csv", sep = ',')
  motif_genes <- motif_genes$HUGO.gene.symbol[motif_genes$HUGO.gene.symbol  %in% rownames(ctstk)]
  apop_genes <- read.csv("gene_lists/apoptosis_hs.csv", sep = ',', header = F)
  apop_genes <- apop_genes$V1[apop_genes$V1  %in% rownames(ctstk)]
  rna_edit_genes <- read.csv("gene_lists/rna_edit_hs.csv", sep = ',', header = F)
  rna_edit_genes <- rna_edit_genes$V1[rna_edit_genes$V1  %in% rownames(ctstk)]
  rna_splicing_genes <- read.csv("gene_lists/rna_splicing_hs.csv", sep = ',', header = F)
  rna_splicing_genes <- rna_splicing_genes$V1[rna_splicing_genes$V1 %in% rownames(ctstk)]
  cytokinesis_genes <- read.csv("gene_lists/cytokinesis.csv", sep = ',', header = F)
  cytokinesis_genes <- cytokinesis_genes$V1[cytokinesis_genes$V1 %in% rownames(ctstk)]
  stemPR_genes <- read.csv("gene_lists/stem_PR.csv", sep = ',', header = F)
  stemPR_genes <- stemPR_genes$V1[stemPR_genes$V1 %in% rownames(ctstk)]
  protein_fold_atp_genes <- read.csv("gene_lists/protein_fold_atp.csv", sep = ',', header = F)
  protein_fold_atp_genes <- protein_fold_atp_genes$V1[protein_fold_atp_genes$V1 %in% rownames(ctstk)]
  dna_dmg_sense_atp_genes <- read.csv("gene_lists/dna_dmg_sense_atp.csv", sep = ',', header = F)
  dna_dmg_sense_atp_genes <- dna_dmg_sense_atp_genes$V1[dna_dmg_sense_atp_genes$V1 %in% rownames(ctstk)]
  element_sensor_gens <- c("SYT1", "EFHB", "EFCAB9", "EFHD1", "SLC39A4") #(calcium , zinc)
  transporter_PR_genes <- read.csv("gene_lists/transporter_PR.csv", sep = ',', header = F)
  transporter_PR_genes <- transporter_PR_genes$V1[transporter_PR_genes$V1 %in% rownames(ctstk)]
  m2m_genes <- c("EZR", "MSN", "ROCK1", "VACM1", "ICAM1", "PEX16", "STXBP3", "PEX26", "RAB7A", "SNX3")
  m2m_genes <- m2m_genes[m2m_genes%in% rownames(ctstk)]
  demethylation_genes <- read.csv("gene_lists/demethylation.csv", sep = ',', header = F)
  demethylation_genes <- demethylation_genes$V1[demethylation_genes$V1 %in% rownames(ctstk)]
  collagen_recptor_genes <- c("OSCAR", "DDR1", "ITGA2", "ITGA11", "DDR2")
  collagen_recptor_genes <- collagen_recptor_genes[ collagen_recptor_genes%in% rownames(ctstk)]
  mlipid_synth_genes <- read.csv("gene_lists/membrane_lipid_synth.csv", sep = ',', header = F)
  mlipid_synth_genes <- unique(mlipid_synth_genes$V1[mlipid_synth_genes$V1 %in% rownames(ctstk)])
  splicosome_genes <- read.csv("gene_lists/splicosome_genes.csv", sep = ';', header = T)
  splicosome_genes <- unique(splicosome_genes$Approved.symbol[splicosome_genes$Approved.symbol %in% rownames(ctstk)])
  splicosome_genes <- c(splicosome_genes, "DGCR14", "CXorf56", "RNU4ATAC", "RNU6ATAC")
  
  hypoxia_genes <- read.csv("gene_lists/hypoxia.csv", sep = ',', header = F)
  hypoxia_genes <- unique(hypoxia_genes$V1[hypoxia_genes$V1 %in% rownames(ctstk)])
  
  att <- as.data.frame(read_tsv('gene_lists/hallmark_hypoxia.tsv',col_types = 'c' ))
  rownames(att) <- att$STANDARD_NAME
  att$STANDARD_NAME <- NULL
  hypoxia_hallmark_genes <- unlist(strsplit(att['GENE_SYMBOLS',], ','))
  hypoxia_hallmark_genes <- hypoxia_hallmark_genes[hypoxia_hallmark_genes %in% rownames(ctstk)]
  
  kidney_genes <- read.csv("gene_lists/kidney genes from literature.csv", sep = ',')
  kidney_genes <- unique(kidney_genes$HUGO.Gene.Symbol[kidney_genes$HUGO.Gene.Symbol %in% rownames(ctstk)])
  
  #splicosome gene lists
  a_complex <- read.csv("gene_lists/a_complex.csv", sep = ';', header = T)
  a_complex <- unique(a_complex$Approved.symbol[a_complex$Approved.symbol %in% rownames(ctstk)])
  b_complex <- read.csv("gene_lists/b_complex.csv", sep = ';', header = T)
  b_complex <- unique(b_complex$Approved.symbol[b_complex$Approved.symbol %in% rownames(ctstk)])
  bact_complex <- read.csv("gene_lists/bact_complex.csv", sep = ';', header = T)
  bact_complex <- unique(bact_complex$Approved.symbol[bact_complex$Approved.symbol %in% rownames(ctstk)])
  c_complex <- read.csv("gene_lists/c_complex.csv", sep = ';', header = T)
  c_complex <- unique(c_complex$Approved.symbol[c_complex$Approved.symbol %in% rownames(ctstk)])
  e_complex <- read.csv("gene_lists/e_complex.csv", sep = ';', header = T)
  e_complex <- unique(e_complex$Approved.symbol[e_complex$Approved.symbol %in% rownames(ctstk)])
  p_complex <- read.csv("gene_lists/p_complex.csv", sep = ';', header = T)
  p_complex <- unique(p_complex$Approved.symbol[p_complex$Approved.symbol %in% rownames(ctstk)])
  minor_complex <- c("RNU4ATAC", "RNU6ATAC")
  
  #genes by expression rank
  down_thres <- -1
  up_thres <- 1
  down_reg_genes <- rownames(wrs_res[order(wrs_res$log2foldChange),])
  down_twice <- which((wrs_res[order(wrs_res$log2foldChange),]$log2foldChange > down_thres) == TRUE)[1]
  up_reg_genes <- rownames(wrs_res[order(wrs_res$log2foldChange, decreasing = T),])
  up_twice <- which((wrs_res[order(wrs_res$log2foldChange, decreasing = T),]$log2foldChange < up_thres) == TRUE)[1]
  
  gene_lists <- list ( "Mesenchymal promoter genes" = mesenchymal_promoter_genes,
                       "Motif genes" = motif_genes,
                       "Apoptosis genes" = apop_genes,
                       "RNA editing genes" = rna_edit_genes,
                       "RNA splicing genes" = rna_splicing_genes,
                       "Cytokinesis genes" = cytokinesis_genes,
                       "Stem cells regulation (+)" = stemPR_genes,
                       "Protein folding Chap' (by ATP)" = protein_fold_atp_genes,
                       "DNA damage sensor (by ATP)"= dna_dmg_sense_atp_genes,
                       "Element sensors" = element_sensor_gens,
                       "Transporter regulation (+)" = transporter_PR_genes,
                       "Membrane/ Protein to Membrane dockning" = m2m_genes,
                       "Demethylation" = demethylation_genes,
                       "Collagen receptor activity" = collagen_recptor_genes,
                       "Membrane lipid biosynthesis" = mlipid_synth_genes,
                       "Hypoxia" = hypoxia_genes,
                       "Hypoxia hallmark" = hypoxia_hallmark_genes,
                       "All" = rownames(ctstk),
                       "Kidney" = kidney_genes,
                       "Splicosome" = splicosome_genes,
                       'Complex A' = a_complex,
                       'Complex B' = b_complex,
                       'Complex Bact' = bact_complex,
                       'Complex C' = c_complex,
                       'Complex E' = e_complex,
                       'Complex P' = p_complex,
                       'Minor Complex' = minor_complex,
                       "Thrombus expression (-)" = down_reg_genes,
                       "Thrombus expression (+)" = up_reg_genes,
                       "Thrombus expression (-/+)" = c(down_reg_genes[1:down_twice], up_reg_genes[1:up_twice]))
  
  gene_lists
}

meta_table_gen <- function(meta_path, meta_more_path)
{
  #load meta data
  meta <- read.csv(meta_path, sep=',')
  meta$Normal <-     ifelse( grepl("N",    meta$Sample.Name), "Yes", "No")
  meta$Thrombosis <- ifelse( grepl("[c,t]",meta$Sample.Name), "Yes", "No")
  meta$State <- "Normal"
  meta$State [meta$Normal == "No" ] <- "Cancer"
  meta$State [grepl("t",meta$Sample.Name)] <- "Thrombus"
  meta$Weird <- ifelse(grepl("R",    meta$Sample.Name), "*", "-")
  
  meta$Thrombosis <- as.factor(meta$Thrombosis)
  meta$State <- as.factor(meta$State)
  meta$Normal <- as.factor(meta$Normal)
  
  meta_more <- read.csv(meta_more_path, sep=',')
  meta_more$T.stage <- as.factor(meta_more$T.stage)
  meta_more$N.stage <- as.factor(meta_more$N.stage)
  meta_more$M.stage <- as.factor(meta_more$M.stage)
  meta_more$Thrombus.level.Mayo. <- as.factor(meta_more$Thrombus.level.Mayo.)
  colnames(meta_more)[colnames(meta_more) == "Thrombus.level.Mayo."] <- "MayoLevel"
  
  meta <- merge(meta, meta_more,all = T)
  
  #order by SRR number
  meta <- meta[order (meta[,"Run"]),]
  
  meta
}

wrs_res_gen <- function(ctstk, meta)
  #wilcoxon rank sum test for differential expressed genes
{
  if(file.exists("data/matrices/wrs_res.csv")){wrs_res <- read.csv(file = "data/matrices/wrs_res.csv", row.names = 1)}else{
    conditions <- as.factor(meta[meta$Normal == "No",]$Thrombosis)
    y <- DGEList(counts=ctstk[, meta$Normal == "No"],group=conditions) #create a group of Yes and No of thrombosis
    ##Remove rows consistently have zero or very low counts
    keep <- filterByExpr(y)
    y <- y[keep,keep.lib.sizes=FALSE]
    ##Perform TMM normalization and transfer to CPM (Counts Per Million)
    y <- calcNormFactors(y,method="TMM")
    count_norm=cpm(y)
    count_norm<-as.data.frame(count_norm)
    
    #Run the Wilcoxon rank-sum test for each gene
    pvalues <- sapply(1:nrow(count_norm),function(i){
      data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
      p=wilcox.test(gene~conditions, data)$p.value
      return(p)
    })
    fdr=p.adjust(pvalues,method = "fdr")
    
    #Calculate the fold-change for each gene
    conditionsLevel<-levels(conditions)
    dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
    dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
    foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
    
    #Output results based on FDR threshold
    outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
    rownames(outRst)=rownames(count_norm)
    outRst=na.omit(outRst)
    fdrThres=0.01
    wrs_res <- as.data.frame (outRst[outRst$FDR<fdrThres,])
    write.csv(wrs_res, file = "data/matrices/wrs_res.csv")
  }
}

ann_colors <- function()
{
  ann_color_list <- list( sex =  c( male = "#4575B4", 
                                    female="#D73027"),
                          Thrombosis = c( No ="#FAEBD7",
                                          Yes="#D73027"),
                          State = c( Normal = "#FAEBD7", 
                                     Cancer = "#E88D7F", 
                                     Thrombus="#D73027"),
                          MayoLevel = c('-' = "#FAEBD7", 
                                        '0' = "#E3BBAC", 
                                        '1' = "#CD8C81", 
                                        '2' = "#B75D55", 
                                        '3' = "#A12E2A", 
                                        '4' = "#8B0000"),
                          T.stage = c('1a'= "#FAEBD7",
                                      '1b'= "#E0DBD2",
                                      '2a'= "#C6CCCD",
                                      '2b'= "#ACBEC8",
                                      '3a'= "#93AFC3",
                                      '3b'= "#79A0BE",
                                      '3c'= "#5F91B9",
                                      '4' = "#4682B4"),
                          N.stage = c('0' = "#FAEBD7", 
                                      '1' = "#A0B6C5",
                                      'x' = "#4682B4"),
                          M.stage = c('0' = "#FAEBD7", 
                                      '1' = "#4682B4"))
  ann_color_list
}

ann_colors_tcga <- function()
{
  ann_color_list <- list( State =  c( 'NT' = "#4575B4",
                                      'TP' ="#D73027"),
                          
                          Stage = c( 'Stage I' = "#FAEBD7",
                                     'Stage II' = "#E88D7F",
                                     'Stage III' = "#D73027",
                                     'Stage IV' = "#A12E2A"))
  ann_color_list
}

kidney_colors <- function()
{
  ann_color_list <- list( Type = c("Uninduced mesenchyme" = "#228B22",
                                   "Cap mesenchyme" = "#328B22",
                                   "EMT" = "#428B22",
                                   
                                   "Podocytes" = "#FA8072",
                                   "Proximal tubules" = "#FA9072",
                                   "Loop of Henle" = "#FAA072",
                                   "Distal tubules" = "#FAB072",
                                   
                                   "Collecting duct" = '#87CEEB',
                                   "Collecting duct (Intercalated cells)" = '#9ACEEB',
                                   "Epithelial" = '#ADCEEB',
                                   
                                   "Macrophages" = '#8B008B',
                                   "Monocytes" = '#AB008B',
                                   "Neutrophils" = '#CB008B',
                                   
                                   "House keeping gene" = 'gold',
                          
                                   "Endothelial cells" = 'pink',
                          
                                   "Cell cycle" = 'azure'))
  ann_color_list
}

Ann_chm <- function()
{
  HeatmapAnnotation( sex = c('male', 'female'),
                     Thrombosis = c('No', 'Yes'),
                     State = c('Normal', 'Cancer', 'Thrombus'),
                     MayoLevel = c('-','0','1','2','3','4'),
                     T.stage = c('1a', '1b' ,'2a', '2b', '3a', '3b', '3c', '4'),
                     N.stage = c('0', '1', 'x'),
                     M.stage = c('0','1'),
    col = list(sex =  c( 'male' = "#4575B4",
                                         'female'="#D73027"),
                               Thrombosis = c( 'No' ="#FAEBD7",
                                               'Yes'="#D73027"),
                               State = c( 'Normal' = "#FAEBD7",
                                          'Cancer' = "#E88D7F",
                                          'Thrombus'="#D73027"),
                               MayoLevel = c('-' = "#FAEBD7",
                                             '0' = "#E3BBAC",
                                             '1' = "#CD8C81",
                                             '2' = "#B75D55",
                                             '3' = "#A12E2A",
                                             '4' = "#8B0000"),
                               T.stage = c('1a'= "#FAEBD7",
                                           '1b'= "#E0DBD2",
                                           '2a'= "#C6CCCD",
                                           '2b'= "#ACBEC8",
                                           '3a'= "#93AFC3",
                                           '3b'= "#79A0BE",
                                           '3c'= "#5F91B9",
                                           '4' = "#4682B4"),
                               N.stage = c('0' = "#FAEBD7",
                                           '1' = "#A0B6C5",
                                           'x' = "#4682B4"),
                               M.stage = c('0' = "#FAEBD7",
                                           '1' = "#4682B4")))
}

ggsashimi_arrow <- function(df_path, dir)
{
  df <- read.csv(df_path, row.names = 1)
  sashimi_list <- list.files(dir, pattern = "*.png")
  df_row_names <- sapply(row.names(df), function(name) {substr(name,0,unlist(gregexpr('_', name))[5]-1)})
  df_row_names <- paste0(df_row_names,'.png')
  for (sashimi in sashimi_list)
  {
    sashimi_image <- image_read(paste0(dir,sashimi))
    sashimi_index <- which (df_row_names == sashimi )[1]
    ex <- (df$exonEnd[sashimi_index] + df$exonStart_0base[sashimi_index])/2
    es<- df$upstreamES[sashimi_index]
    ee<- df$downstreamEE[sashimi_index]
    arw_offset <- image_info(sashimi_image)$width * (0.983 + 0.1235*(ee-ex)/(ex-es))/(1 + (ee-ex)/(ex-es))
    pos_string <- paste0('+',arw_offset,'+60')
    # 0.1235 - x1%, 0.985 - x2%
    sashimi_image <- image_annotate(sashimi_image,sprintf('\u2193'), 
                                    size = 60, location = pos_string,color = "forestgreen")
    ggsave(paste0(dir,sashimi),
           rasterGrob(as.raster(sashimi_image)), 
           width = image_info(sashimi_image)$width, 
           height = image_info(sashimi_image)$height, 
           units = "px")
  }
}

feature_plot <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                         point_size = 7, #point size
                         min_point_size = 0.01,
                         text_size = 15,#size of label text
                         inclusion_level = 1, #for feature plot need vector of incLevels
                         myTitle = "Graph",
                         color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                         name_t = "Inclusion level",
                         log2 = FALSE)
{
  if (log2 == TRUE)
  {
    inclusion_level <- log2 (inclusion_level + 1)
  }
 # else {
    inclusion_level <- (inclusion_level - min(inclusion_level)) / (max(inclusion_level) - min(inclusion_level))
    #}
  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2])
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")
  
  color_breaks <- quantile_breaks(inclusion_level,101)
  point_style <- geom_point(aes(fill = inclusion_level), shape = 21, colour = 'black', 
                            size=point_size * (inclusion_level+min_point_size))
                            #size=point_size * log(100*(inclusion_level+0.0101)))
  
  gradient_style <- scale_fill_gradientn(colors = colorRampPalette(color_wheel)(length(color_breaks)-1),
                                        name = name_t, values = color_breaks )
  
  ggplot(data=d, aes(x=d$PC1, y=d$PC2)) +
    point_style + gradient_style + xlab(xs) + ylab(ys) +
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}

pca_plot <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                     point_size = 7, #point size
                     text_size = 15,#size of label text
                     myTitle = "Graph",
                     groupMeta = NULL,
                     exFactor = NULL,
                     meta_type = 0,
                     color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                     markA = NULL,
                     markB = NULL)
{
  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  if (!is.null(groupMeta))
  {
    if (meta_type == 0)
    group <- factor(patient_state(groupMeta), 
                   levels = c('Normal', 
                              'Cancer', 
                              'Primary Tumor', 
                              'Thrombus Tumor'))
    else
    {
      group1 <- factor(patient_state(groupMeta[1:189,]), 
                      levels = c('Normal', 
                                 'Cancer', 
                                 'Primary Tumor', 
                                 'Thrombus Tumor'))
      group2 <- factor(exFactor)
      group <- factor (c(group1,group2))
    }
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2], State = group)
  
  if(is.null(markA)) markA <- geom_blank()
  else{markA <- geom_point(data = d[markA,],aes(x = PC1, y= PC2),pch=21,
                              fill=NA, size=point_size+4, colour="#00C08B", stroke=2)}
  
  if(is.null(markB)) markB <- geom_blank()
  else {markB <- geom_point(data = d[markB,],aes(x = PC1, y= PC2),pch=21, 
                               fill=NA, size=point_size+4, colour="#FF64B0", stroke=2)}
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")

  ggplot(data=d, aes(x=PC1, y=PC2)) +
           geom_point(aes(fill = d$State), shape = 21, colour = 'black', size=point_size) + 
    markA + markB + xlab(xs) + ylab(ys) +
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}


factor_plot <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                     point_size = 7, #point size
                     text_size = 15,#size of label text
                     myTitle = "Graph",
                     groupMeta = NULL,
                     color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                     markA = NULL,
                     markB = NULL,
                     fill_title = 'Levels')
{
  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  if (!is.null(groupMeta))
    group <-factor(groupMeta)
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2], State = group)
  
  if(is.null(markA)) markA <- geom_blank()
  else{markA <- geom_point(data = d[markA,],aes(x = PC1, y= PC2),pch=21,
                           fill=NA, size=point_size+4, colour="#00C08B", stroke=2)}
  
  if(is.null(markB)) markB <- geom_blank()
  else {markB <- geom_point(data = d[markB,],aes(x = PC1, y= PC2),pch=21, 
                            fill=NA, size=point_size+4, colour="#FF64B0", stroke=2)}
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")
  
  ggplot(data=d, aes(x=PC1, y=PC2)) +
    geom_point(aes(fill = d$State), shape = 21, colour = 'black', size=point_size) + 
    markA + markB + xlab(xs) + ylab(ys) + labs(fill = fill_title) + 
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}

age_plot <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                         point_size = 7, #point size
                         min_point_size = 0.01,
                         text_size = 15,#size of label text
                         inclusion_level_o = 1, #for feature plot need vector of incLevels
                         myTitle = "Graph",
                         color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                         name_t = "Inclusion level")
{
  inclusion_level <- (inclusion_level_o - min(inclusion_level_o)) / (max(inclusion_level_o) - min(inclusion_level_o))
  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2])
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")
  
  color_breaks <- quantile_breaks(inclusion_level,101)
  point_style <- geom_point(aes(fill = inclusion_level), shape = 21, colour = 'black', 
                            size=point_size * (inclusion_level+min_point_size))
  
   gradient_style <- scale_fill_gradientn(colors = colorRampPalette(color_wheel)(length(color_breaks)-1),
                                          name = name_t, 
                                          values = color_breaks, 
                                          breaks = c(min(inclusion_level), median(inclusion_level) ,max(inclusion_level)),
                                          labels=  c(min(inclusion_level_o), median(inclusion_level_o) ,max(inclusion_level_o)))
  
  ggplot(data=d, aes(x=d$PC1, y=d$PC2)) +
    point_style  + gradient_style + 
    xlab(xs) + ylab(ys) +
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}





feature_plot2 <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                         point_size = 7, #point size
                         min_point_size = 0.01,
                         text_size = 15,#size of label text
                         inclusion_level = 1, #for feature plot need vector of incLevels
                         myTitle = "Graph",
                         color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                         name_t = "Inclusion level",
                         log2 = FALSE)
{
  if (log2 == TRUE)
  {
    inclusion_level <- log2 (inclusion_level + 1)
  }

  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2])
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")
  
  color_breaks <- quantile_breaks2(inclusion_level)#seq(min(inclusion_level), max(inclusion_level), length.out = 10)
  point_style <- geom_point(aes(fill = inclusion_level), shape = 21, colour = 'black', 
                            size=point_size * ((inclusion_level - min(inclusion_level)) / (max(inclusion_level) - min(inclusion_level))+min_point_size))
  #size=point_size * log(100*(inclusion_level+0.0101)))
  
  gradient_style <- scale_fill_gradientn(colors = colorRampPalette(color_wheel)(length(color_breaks)-1),
                                         name = name_t, 
                                         #breaks = c(min(color_breaks), max(color_breaks)),
                                         breaks = color_breaks,
                                         labels = formatC(c(min(color_breaks), max(color_breaks)), digits = 2, format = 'f'))
  
  ggplot(data=d, aes(x=d$PC1, y=d$PC2)) +
    point_style + gradient_style + xlab(xs) + ylab(ys) +
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}

quantile_breaks3 <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = TRUE)
  breaks[!duplicated(breaks)]
}

feature_plot3 <- function(object,pc1 = 1, pc2 = 2, #pcs to use
                     point_size = 7, #point size
                     min_point_size = 0.01,
                     text_size = 15,#size of label text
                     inclusion_level_o = 1, #for feature plot need vector of incLevels
                     myTitle = "Graph",
                     color_wheel = c("#4575B4","#FFFFBF", "#D73027"),
                     name_t = "Inclusion level")
{
  inclusion_level <- (inclusion_level_o - min(inclusion_level_o)) / (max(inclusion_level_o) - min(inclusion_level_o))

  object <- object[apply(object, 1, var, na.rm=TRUE) != 0,]
  df_pca <- t(object)
  pca <-  prcomp(df_pca[,apply(df_pca,2,function(x) sum(x)!=0)], scale. = T, center = T)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )   # the contribution to the total variance for each component
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2])
  
  #create pc labels
  xs = paste0("PC", as.character(pc1),": ", round(percentVar[pc1] * 100),"% variance")
  ys = paste0("PC", as.character(pc2),": ", round(percentVar[pc2] * 100),"% variance")
  
  color_breaks <- quantile_breaks(inclusion_level,101) #inclusion_level
  point_style <- geom_point(aes(fill = inclusion_level), shape = 21, colour = 'black', 
                            size=point_size * (inclusion_level+min_point_size))
  
  gradient_style <- scale_fill_gradientn(colors = colorRampPalette(color_wheel)(length(color_breaks)-1),
                                         name = name_t, 
                                         #values = color_breaks, 
                                         #breaks = c(min(inclusion_level),max(inclusion_level)),
                                         #labels=  formatC(c(min(inclusion_level_o),max(inclusion_level_o)),
                                        #                  digits = 2, format = 'f')
                                        # )
                                         values = color_breaks, 
                                         breaks = c(min(inclusion_level), median(inclusion_level) ,max(inclusion_level)),
                                         labels=  formatC2(c(min(inclusion_level_o), median(inclusion_level_o) ,max(inclusion_level_o))))#, digits = 2, format = 'f'))
  
  
  
  ggplot(data=d, aes(x=d$PC1, y=d$PC2)) +
    point_style  + gradient_style + 
    xlab(xs) + ylab(ys) +
    coord_fixed() + ggtitle(myTitle) + 
    theme(plot.title = element_text(size=text_size),
          legend.title=element_text(size=text_size))
}


formatC2 <- function(inclusion_level_o)
{
  x <- formatC(c(min(inclusion_level_o), median(inclusion_level_o) ,max(inclusion_level_o)), digits = 2, format = 'f')
  sapply(x, function(i) if (i == '0.00') '0' else i)
}
