#---------------------------------------------------------------#
# Axis clusters                                                 #
# Figure 3 code                                                 #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
# ggplot settings
ggplot2_theme <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=10),
        plot.subtitle=element_text(hjust=0.5, size=10),
        axis.text=element_text(colour='black'),
        axis.ticks=element_line(colour='black'))
theme_set(ggplot2_theme)

#----------------------------------------------------------------#
# Figure 3 A-B                                                   #
#----------------------------------------------------------------#
# function for average signal in regions
getmean <- function(regions, signalfile) {
  regionsignals = vector()
  for(i in 1:length(regions)){
    regionsignals[i] = regioneR::meanInRegions(A=regions[i], x=signalfile)
  }
  return(regionsignals)
}
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

# A: H3
h3 = import_bedGraph('H3_Zhu2015-SK1Yue-PM_B3-MACS2_FE.bdg')
genAvg <- hwglabr2::average_chr_signal(h3)$genome_avrg
h3$score <- h3$score/genAvg

h3_in_cluster = getmean(clusters,h3)
h3_in_desert = getmean(deserts,h3)
hist(h3_in_cluster)
hist(h3_in_desert)

h3_in_cluster_df = data.frame(region = 'cluster',signal = h3_in_cluster)
h3_in_desert_df = data.frame(region = 'desert',signal = h3_in_desert)
h3_all = rbind(h3_in_cluster_df,h3_in_desert_df)
hist(h3_all$signal)

p <- ggplot(h3_all, aes(x=region, y=signal, fill = region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") +
  ylab("H3 signal") + xlab("") + geom_boxplot(fill=NA)
p
wilcox.test(h3_in_cluster,h3_in_desert) #p-value < 2.2e-16

# A: H4
h4 = import_bedGraph('H4_Hu2015_4h-SK1Yue-PM_B3-MACS2_FE.bdg')
genAvg <- hwglabr2::average_chr_signal(h4)$genome_avrg
h4$score <- h4$score/genAvg

h4_in_cluster = getmean(clusters,h4)
h4_in_desert = getmean(deserts,h4)
hist(h4_in_cluster)
hist(h4_in_desert)

h4_in_cluster_df = data.frame(region = 'cluster',signal = h4_in_cluster)
h4_in_desert_df = data.frame(region = 'desert',signal = h4_in_desert)
h4_all = rbind(h4_in_cluster_df,h4_in_desert_df)
hist(h4_all$signal)

p <- ggplot(h4_all, aes(x=region, y=signal, fill = region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") +
  ylab("H4 signal") + xlab("") + geom_boxplot(fill=NA)
p
wilcox.test(h4_in_cluster,h4_in_desert) #p-value = 0.0003795

#----------------------------------------------------------------#
# Figure 3 C                                                     #
#----------------------------------------------------------------#
mnase = hwglabr2::import_bedGraph('Nucleosome_reps-SK1-MACS2_treat_pileup.bdg')
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
mnase = gendiv(mnase)

clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
gff = get_gff("SK1Yue")
gff = gff[gff$type == 'gene']
# find genes that overlap with cluters or deserts
hits <- findOverlaps(query = clusters,subject = gff)
overlaps <- pintersect(clusters[queryHits(hits)], gff[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
cluster_gff <- gff[subjectHits(hits)]
rm(hits);rm(overlaps)

hits <- findOverlaps(deserts,gff)
overlaps <- pintersect(deserts[queryHits(hits)], gff[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
deserts_gff <- gff[subjectHits(hits)]

# desert genes
deserts_gffplus = deserts_gff[strand(deserts_gff)=='+']
end(deserts_gffplus) = start(deserts_gffplus)
desert_gffminus = deserts_gff[strand(deserts_gff)=='-']
start(desert_gffminus) = end(desert_gffminus)
deserts_gffall = c(deserts_gffplus,desert_gffminus)

Mnase_desertgffall = EnrichedHeatmap::normalizeToMatrix(mnase,
                                                        deserts_gffall,
                                                        extend=c(100,2000), w=1,
                                                        mean_mode="weighted",
                                                        value_column="score")
Mnase_desertgffall_avrg <- hwglabr2::signal_mean_and_ci(Mnase_desertgffall,
                                                        ci=0.95, rep_bootstrap=1000,
                                                        na_rm=TRUE)
#cluster genes
cluster_gffplus = cluster_gff[strand(cluster_gff)=='+']
end(cluster_gffplus) = start(cluster_gffplus)
cluster_gffminus = cluster_gff[strand(cluster_gff)=='-']
start(cluster_gffminus) = end(cluster_gffminus)
cluster_gffall = c(cluster_gffplus,cluster_gffminus)

Mnase_clustergffall = EnrichedHeatmap::normalizeToMatrix(mnase,
                                                         cluster_gffall,
                                                         extend=c(100,2000), w=1,
                                                         mean_mode="weighted",
                                                         value_column="score")
Mnase_clustergffall_avrg <- hwglabr2::signal_mean_and_ci(Mnase_clustergffplus,
                                                         ci=0.95, rep_bootstrap=1000,
                                                         na_rm=TRUE)

Mnase_desertgffall_avrgdf <- data.frame(Position=seq(-100, 2000), Region = 'desert',Mnase_desertgffall_avrg)
Mnase_clustergffall_avrgdf <- data.frame(Position=seq(-100, 2000), Region = 'cluster',Mnase_clustergffall_avrg)
allstrand = rbind(Mnase_desertgffall_avrgdf,Mnase_clustergffall_avrgdf)

p <- ggplot(allstrand, aes(Position, Mean, fill = Region, color = Region)) +
  labs(title = "Signal at ATGs",
       x = "Distance to ATG (bp)", y = "Average MNase\nsignal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-100, 0, 2000),
                     labels = c("-100bp", "ATG", "2000bp"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA)
p + geom_line() + ylim(0.8,1.55)

#----------------------------------------------------------------#
# Figure 3D                                                      #
#----------------------------------------------------------------#
clusters = rtracklayer::import.bed('\clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
hits = findOverlaps(clusters,hotspots)
overlaps <- pintersect(clusters[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_cluster <- hotspots[subjectHits(hits)]
rm(hits);rm(overlaps)

hits = findOverlaps(deserts,hotspots)
overlaps <- pintersect(deserts[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_desert <- hotspots[subjectHits(hits)]

mcols(hotspots_desert)['class'] = 'desert'
mcols(hotspots_cluster)['class'] = 'cluster'

hotspots_desert_sort = hotspots_desert[order(width(hotspots_desert),decreasing = T)]
hotspots_cluster_sort = hotspots_cluster[order(width(hotspots_cluster),decreasing = T)]
hotspots_all = c(hotspots_desert_sort,hotspots_cluster_sort)
midpoint <- floor(width(hotspots_all) / 2)
start(hotspots_all) <- start(hotspots_all) + midpoint
end(hotspots_all) <- start(hotspots_all)

if(!(dir.exists("hotspots_pdf"))) {
  dir.create("hotspots_pdf")
}

# bedgraph files to analyze
# File in folder (clusterdesert)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  Red1_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Red1_bg)$genome_avrg
  Red1_bg$score <- Red1_bg$score/genAvg

  mat1 <- normalizeToMatrix(Red1_bg, hotspots_all, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  col_fun <- circlize::colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))

  pdf(paste0("hotspots_pdf/Nucleosomes_around_hotspots_heatmap_sort.pdf"))
  print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c('green','red')),
                                                                                 show_error = TRUE,pos_line=FALSE)),
                        row_title_rot = 0,#	top_annotation_height = unit(2, "cm"),
                        axis_name = c("-1 kb", "hotspots", "1 kb"),
                        row_order = 1:length(hotspots_all),
                        split=hotspots_all$class,
                        column_title ='Nucleosomes'))
  dev.off()
}


#----------------------------------------------------------------#
# Figure 3D                                                      #
#----------------------------------------------------------------#
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")

hits = findOverlaps(clusters,hotspots)
overlaps <- pintersect(clusters[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_cluster <- hotspots[subjectHits(hits)]
rm(hits);rm(overlaps)

hits = findOverlaps(deserts,hotspots)
overlaps <- pintersect(deserts[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_desert <- hotspots[subjectHits(hits)]

mcols(hotspots_desert)['class'] = 'desert'
mcols(hotspots_cluster)['class'] = 'cluster'

hotspots_desert_sort = hotspots_desert[order(width(hotspots_desert),decreasing = T)]
hotspots_cluster_sort = hotspots_cluster[order(width(hotspots_cluster),decreasing = T)]
hotspots_all = c(hotspots_desert_sort,hotspots_cluster_sort)
midpoint <- floor(width(hotspots_all) / 2)
start(hotspots_all) <- start(hotspots_all) + midpoint
end(hotspots_all) <- start(hotspots_all)

dir.create("hotspots_pdf")
bedgraphs <- "Nucleosome_reps-SK1-MACS2_treat_pileup.bdg"

# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  Red1_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Red1_bg)$genome_avrg
  Red1_bg$score <- Red1_bg$score/genAvg

  mat1 <- normalizeToMatrix(Red1_bg, hotspots_all, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  col_fun <- circlize::colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))

  pdf(paste0("hotspots_pdf/Nucleosomes_reps_around_hotspots_heatmap_sort.pdf"))
  print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c('green','red')),
                                                                                 show_error = TRUE,pos_line=FALSE)),
                        row_title_rot = 0,
                        axis_name = c("-1 kb", "hotspots", "1 kb"),
                        row_order = 1:length(hotspots_all),
                        split=hotspots_all$class,
                        column_title ="Nucleosomes"))
  dev.off()
}

#----------------------------------------------------------------#
# Figure 3 E-H                                                   #
#----------------------------------------------------------------#
# Divide hotspots into cluster or desert
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
hits = findOverlaps(clusters,hotspots)
overlaps <- pintersect(clusters[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_cluster <- hotspots[subjectHits(hits)]
rm(hits);rm(overlaps)

hits = findOverlaps(deserts,hotspots)
overlaps <- pintersect(deserts[queryHits(hits)], hotspots[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(hotspots[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
hotspots_desert <- hotspots[subjectHits(hits)]

mcols(hotspots_desert)['class'] = 'desert'
mcols(hotspots_cluster)['class'] = 'cluster'

hotspots_all = c(hotspots_desert,hotspots_cluster)
midpoint <- floor(width(hotspots_all) / 2)
start(hotspots_all) <- start(hotspots_all) + midpoint
end(hotspots_all) <- start(hotspots_all)

# Put the following files into folder called "clusterdesert_hotspots"
# E: H4K44ac_Hu2015_4h-SK1Yue-PM_B3-MACS2_FE.bdg
# F: H4_Hu2015_4h-SK1Yue-PM_B3-MACS2_FE.bdg
# G: H3K4me3_Zhu2015-SK1Yue-PM_B3-MACS2_FE.bdg
# H: H3_Zhu2015-SK1Yue-PM_B3-MACS2_FE.bdg

bedgraphs <- list.files(path="~/Desktop/clusterdesert_hotspots",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

dir.create("hotspots_pdf")
# function to plot average signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][6],split="_")[[1]][1]

  signalfile <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(signalfile)$genome_avrg
  signalfile$score <- signalfile$score/genAvg

  mat1 <- normalizeToMatrix(signalfile, hotspots_all[hotspots_all$class=='cluster'], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)

  mat2 <- normalizeToMatrix(signalfile, hotspots_all[hotspots_all$class=='desert'], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)

  mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
  mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
  alldata = rbind(mat1_avrg_df,mat2_avrg_df)

  library(ggplot2)
  pdf(paste0("hotspots_pdf/",naming,"_around_hotspots_avg.pdf"), width = 6, height = 4)
  p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
    labs(title = "Signal at hotspots",
         x = "Distance to hotspot (bp)", y = "Average\nChIP-seq signal") +
    geom_vline(xintercept = 0, lty = 3) +
    scale_x_continuous(breaks = c(-199, 0, 200),
                       labels = c("-1 kb", "hotspots", "1 kb"))
  p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)
  print(p + geom_line())

  dev.off()
}
