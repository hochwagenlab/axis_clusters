#---------------------------------------------------------------#
# Axis Islands -- Islands = clusters in code                    #
# Figure 3 code                                                 #
#---------------------------------------------------------------#

library(GenomicRanges)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
library(circlize)

#Figure 3A
setwd('/Users/darmokandjalad/Documents/HI-Scripts_Analysis/IslandPaper/GitHub')

clusters_sort <- clusters[order(width(clusters),decreasing = T)]

flank <- floor(GenomicRanges::width(clusters_sort)/2)
GenomicRanges::start(clusters_sort) <- GenomicRanges::start(clusters_sort) - flank
GenomicRanges::end(clusters_sort) <- GenomicRanges::end(clusters_sort) + flank
clusters_sort <- GenomicRanges::trim(clusters_sort)
mcols(clusters_sort) <- DataFrame(class=c(rep(1:3, each=length(clusters_sort)/3)))

Red1_rec8 <- hwglabr2::import_bedGraph("Red1-rec8D-39-62-193-90-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_rec8d = gendiv(Red1_rec8)


# Concert window_size to n of windows
n_windows <- floor(1000 / 5)
Red1_clustermat <- EnrichedHeatmap::normalizeToMatrix(Red1_rec8d, clusters_sort,
                                                      value_column="score",
                                                      mean_mode="weighted",
                                                      extend=0, k=n_windows,
                                                      empty_value=NA)


col_fun <- colorRamp2(quantile(Red1_clustermat, c(0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("gold1", "lightgoldenrod1","white", "lightsteelblue1","purple4"))
partition <- clusters_sort$class
EnrichedHeatmap(Red1_clustermat, col = col_fun, name = "Red1 in rec8", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:3),show_error = TRUE)),
                row_order = 1:length(clusters_sort),
                split=clusters_sort$class,pos_line = T,
                axis_name = c("","5' ","3' ",""))+
  Heatmap(partition, col = structure(1:3, names = as.character(1:3)), name = " ",row_order = 1:length(clusters_sort),
          show_row_names = FALSE, width = unit(5, "mm"))


#----------------------------------------------------------------#
# Figure 3B,C,D                                                  #
#----------------------------------------------------------------#
# Divide axis sites into cluster or desert
setwd('/Users/darmokandjalad/Desktop/clusterdesert/axis')
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



# B: ssc2 AH8867B-379-767-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
#    smc4 AH6408I-144-183-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# C: top1-13myc AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz
#    Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# D: Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak
#    9120-antiRed1-20170817-20171115-20171207-20180530-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz


bedgraphs <- list.files(path="~/Desktop/clusterdesert/axis",pattern="_FE.bdg",recursive=T,full.names=T)
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
