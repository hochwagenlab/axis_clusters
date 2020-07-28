#---------------------------------------------------------------#
# Axis clusters                                                 #
# Supplemental Figure 1 code                                    #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
#---------------------------------------------------------------#
# Fig S1 A-B                                                    #
#---------------------------------------------------------------#
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

# Fig S1 A
Red1_WT = import_bedGraph("Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_WTd = gendiv(Red1_WT)

signal_desertgff <- hwglabr2::signal_at_orf2(signal_data=Red1_WTd, gff=deserts_gff,
                                             write_to_file=FALSE)
signal_desertgff_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=signal_desertgff,
                                                         ci=0.95, rep_bootstrap=1000,
                                                         na_rm=TRUE)
signal_clustergff <- hwglabr2::signal_at_orf2(signal_data=Red1_WTd, gff=cluster_gff,
                                              write_to_file=FALSE)
signal_clustergff_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=signal_clustergff,
                                                          ci=0.95, rep_bootstrap=1000,
                                                          na_rm=TRUE)
signal_desertgff_metaORFdf <- data.frame(Region='desert',Position=seq(1, 1000), signal_desertgff_metaORF)
signal_clustergff_metaORFdf <- data.frame(Region='cluster',Position=seq(1, 1000), signal_clustergff_metaORF)
allwtgroups <- rbind(signal_desertgff_metaORFdf,signal_clustergff_metaORFdf)

# Set up the plot
library(ggplot2)
p <- ggplot(title='wt',allwtgroups, aes(Position, Mean, group = Region,fill=Region,colour=Region)) +
  labs(x = "Scaled ORF", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = c(250, 750), lty = 3) +
  scale_x_continuous(breaks = c(0, 250, 750, 1000),
                     labels = c("", "start", "stop", ""))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,color=NA) +
  geom_line()
p

#--------------------------------------------------------------#
# Fig S1 B
Red1_rec8 <- hwglabr2::import_bedGraph("Red1-rec8D-39-62-193-90-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_rec8d = gendiv(Red1_rec8)

Red1rec8_desertgff <- hwglabr2::signal_at_orf2(signal_data=Red1_rec8d, gff=deserts_gff,
                                               write_to_file=FALSE)
Red1rec8_desertgff_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Red1rec8_desertgff,
                                                           ci=0.95, rep_bootstrap=1000,
                                                           na_rm=TRUE)
Red1rec8_clustergff <- hwglabr2::signal_at_orf2(signal_data=Red1_rec8d, gff=cluster_gff,
                                                write_to_file=FALSE)
Red1rec8_clustergff_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Red1rec8_clustergff,
                                                            ci=0.95, rep_bootstrap=1000,
                                                            na_rm=TRUE)
Red1rec8_desertgff_metaORFdf <- data.frame(Region='desert',Position=seq(1, 1000), Red1rec8_desertgff_metaORF)
Red1rec8_clustergff_metaORFdf <- data.frame(Region='cluster',Position=seq(1, 1000), Red1rec8_clustergff_metaORF)
allrec8groups <- rbind(Red1rec8_desertgff_metaORFdf,Red1rec8_clustergff_metaORFdf)

# Set up the plot
library(ggplot2)
p <- ggplot(title='rec8',allrec8groups, aes(Position, Mean, group= Region,fill=Region,colour=Region)) +
  labs(x = "Scaled ORF", y = "Red1 rec8\nChIP-seq signal") +
  geom_vline(xintercept = c(250, 750), lty = 3) +
  scale_x_continuous(breaks = c(0, 250, 750, 1000),
                     labels = c("", "start", "stop", ""))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,color=NA) +
  geom_line()
p

#---------------------------------------------------------------#
# Fig S1 C-G                                                    #
#---------------------------------------------------------------#
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
# C: AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz
# D: Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# E: Top2-wildtype-0h-412-503-530-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# F: Top2-spo11YF-514-543-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# G: Top2-spo11D-420-513-541-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz

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

#---------------------------------------------------------------#
# Fig S1 H-K                                                    #
#---------------------------------------------------------------#
# Divide axis sites into cluster or desert
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
axis = hwglabr2::get_Red1_summits("SK1Yue")

hits = findOverlaps(clusters,axis)
axis_cluster <- axis[subjectHits(hits)]
rm(hits)
hits = findOverlaps(deserts,axis)
axis_desert <- axis[subjectHits(hits)]
rm(hits)
subset(axis_cluster, (name %in% axis_desert$name))
mcols(axis_desert)['class'] = 'desert'
mcols(axis_cluster)['class'] = 'cluster'

# Put the following files into folder called "clusterdesert_axis"
# H: Red1-WT-34C-410-495-528-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# I: Red1-top2-4-411-498-535-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# J: Hop1-WT-34C-494-532-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# K: Hop1-top2-4-34C-497-536-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz

bedgraphs <- list.files(path="~/Desktop/clusterdesert/axis",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

dir.create("axis_pdf")
# function to plot average signal at cluster and desert axis sites
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][6],split="_")[[1]][1]

  signalfile <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(signalfile)$genome_avrg
  signalfile$score <- signalfile$score/genAvg

  mat1 <- normalizeToMatrix(signalfile, axis_cluster, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                    ci=0.95, rep_bootstrap=1000,
                    na_rm=TRUE)

  mat2 <- normalizeToMatrix(signalfile, axis_desert, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)

  mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
  mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
  alldata = rbind(mat1_avrg_df,mat2_avrg_df)

  library(ggplot2)
  pdf(paste0("axis_pdf/",naming,"_around_axis_avg.pdf"), width = 6, height = 4)
  p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
    labs(title = "Signal at Red1 summits",
         x = "Distance to axis (bp)", y = "Average\nChIP-seq signal") +
    geom_vline(xintercept = 0, lty = 3) +
    scale_x_continuous(breaks = c(-199, 0, 200),
                       labels = c("-1 kb", "summit", "1 kb"))

  # Add confidence interval as a ribbon
  p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)

  # Add signal line
  print(p + geom_line())

  dev.off()
}
