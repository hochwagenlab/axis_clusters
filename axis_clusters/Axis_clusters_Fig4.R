#---------------------------------------------------------------#
# Axis clusters                                                 #
# Figure 4 code                                                 #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)

#----------------------------------------------------------------#
# Figure 4 A                                                     #
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
bedgraphs <- "Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg"

# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  Red1_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Red1_bg)$genome_avrg
  Red1_bg$score <- Red1_bg$score/genAvg

  mat1 <- normalizeToMatrix(Red1_bg, hotspots_all, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  col_fun <- circlize::colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))

  pdf(paste0("hotspots_pdf/Spo11oligos_around_hotspots_heatmap_sort.pdf"))
  print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c('green','red')),
                                                                                 show_error = TRUE,pos_line=FALSE)),
                        row_title_rot = 0,
                        axis_name = c("-1 kb", "hotspots", "1 kb"),
                        row_order = 1:length(hotspots_all),
                        split=hotspots_all$class,
                        column_title ="Spo11 oligos"))
  dev.off()
}

#----------------------------------------------------------------#
# Figure 4 B,D,F                                                 #
#----------------------------------------------------------------#
clusters = rtracklayer::import.bed('clusters_joined.bed')

# B
spo11oligos = rtracklayer::import.bedGraph('Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg')
genAvg <- hwglabr2::average_chr_signal(spo11oligos)$genome_avrg
spo11oligos$score <- spo11oligos$score/genAvg

clusters_8 = plot_signal(spo11oligos,'VIII',1000)
clusters_8 = data.frame(clusters_8)
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
cen_mid = genome_info[seqnames(genome_info)=="chrVIII"]
cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
max_score = max(spo11oligos$score)
a <- ggplot(clusters_8,aes(position,signal)) +
  geom_line(position='identity') + ylab("chrVIII")
a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                    size = 4, colour = 'green',shape=20)
clusterregion = data.frame(clusters[seqnames(clusters)=="chrVIII"])
a <- a + geom_segment(clusterregion, size = 2,alpha = 0.6,
                      mapping=aes(x = start, y = max_score, xend = end,
                                  yend = max_score, colour = "segment"))
a

# D
pH2A = import_bedGraph('AH6179-pH2A-T3-344-348-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
genAvg <- hwglabr2::average_chr_signal(pH2A)$genome_avrg
pH2A$score <- pH2A$score/genAvg

clusters_8 = plot_signal(pH2A,'VIII',1000)
clusters_8 = data.frame(clusters_8)
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
cen_mid = genome_info[seqnames(genome_info)=="chrVIII"]
cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
max_score = max(pH2A$score)
a <- ggplot(clusters_8,aes(position,signal)) +
  geom_line(position='identity') + ylab("chrVIII")
a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                    size = 4, colour = 'green',shape=20)
clusterregion = data.frame(clusters[seqnames(clusters)=="chrVIII"])
a <- a + geom_segment(clusterregion, size = 2,alpha = 0.6,
                      mapping=aes(x = start, y = max_score, xend = end,
                                  yend = max_score, colour = "segment"))
a

# F
Zip3 = import_bedGraph('Zip3-Flag_DeMuyt2018-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
genAvg <- hwglabr2::average_chr_signal(Zip3)$genome_avrg
Zip3$score <- Zip3$score/genAvg

clusters_8 = plot_signal(Zip3,'VIII',1000)
clusters_8 = data.frame(clusters_8)
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
cen_mid = genome_info[seqnames(genome_info)=="chrVIII"]
cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
max_score = max(Zip3$score)
a <- ggplot(clusters_8,aes(position,signal)) +
  geom_line(position='identity') + ylab("chrVIII")
a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                    size = 4, colour = 'green',shape=20)
clusterregion = data.frame(clusters[seqnames(clusters)=="chrVIII"])
a <- a + geom_segment(clusterregion, size = 2,alpha = 0.6,
                      mapping=aes(x = start, y = max_score, xend = end,
                                  yend = max_score, colour = "segment"))
a


#----------------------------------------------------------------#
# Figure 4 C,E,G                                                 #
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
# C: Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg
# E: AH6179-pH2A-T3-344-348-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# G: Zip3-Flag_DeMuyt2018-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz

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
