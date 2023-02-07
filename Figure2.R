#---------------------------------------------------------------#
# Axis Islands -- Islands = clusters in code                    #
# Figure 2 code                                                 #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
library(patchwork)

# Create working folder with necessary files 
setwd('/Users/darmokandjalad/Documents/HI-Scripts_Analysis/IslandPaper/GitHub')

#----------------------------------------------------------------#
# Figure 2A                                                      #
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
# Figure 2 B,C                                                   #
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


# Figure 2B
bedgraph_file <- 'Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg'
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
p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
  labs(title = "Signal at hotspots",
       x = "Distance to hotspot (bp)", y = "Average\nChIP-seq signal") +
  geom_line(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-1 kb", "hotspots", "1 kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)


#Figure 2C
# G: Zip3-Flag_DeMuyt2018-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
bedgraph_file <- 'Zip3Flag1_trmd_NgNt-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz'

signalfile <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
genAvg <- hwglabr2::average_chr_signal(signalfile)$genome_avrg
signalfile$score <- signalfile$score/genAvg
signalfile <- SMsigN
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
p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
  labs(title = "Signal at hotspots",
       x = "Distance to hotspot (bp)", y = "Average\nChIP-seq signal") +
  geom_line(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-1 kb", "hotspots", "1 kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)