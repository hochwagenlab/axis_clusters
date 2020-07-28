#---------------------------------------------------------------#
# Axis clusters                                                 #
# Figure 1 code                                                 #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
#----------------------------------------------------------------#
# Defining clusters
#----------------------------------------------------------------#
#                                                                #
#          Define clusters with joining nearby regions           #
#                                                                #
#----------------------------------------------------------------#
Red1_rec8 <- hwglabr2::import_bedGraph("Red1-rec8D-39-62-193-90-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_rec8d = gendiv(Red1_rec8)

tilesize=5000
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
rec8ex <- sort(GenomeInfoDb::sortSeqlevels(Red1_rec8d))
GenomeInfoDb::seqlengths(rec8ex) <- GenomeInfoDb::seqlengths(genome_info)
bins_rec8ex <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(rec8ex),
                                         tilewidth=tilesize,
                                         cut.last.tile.in.chrom=TRUE)
score_rec8ex <- GenomicRanges::coverage(rec8ex, weight="score")
bins_rec8ex <- GenomicRanges::binnedAverage(bins_rec8ex, score_rec8ex, "binned_score")

clusters <- bins_rec8ex[bins_rec8ex$binned_score>=1.75* sd(bins_rec8ex$binned_score)]
deserts <- bins_rec8ex[bins_rec8ex$binned_score<1.75* sd(bins_rec8ex$binned_score)]

joinclusters = regioneR::joinRegions(A = clusters, min.dist = 200) #join regions 200 bp or less apart
joindeserts = regioneR::joinRegions(A = deserts, min.dist = 200) #join regions 200 bp or less apart

write.table(data.frame(granges(joinclusters))[,1:3],'clusters_joined.bed',quote = F,row.names = F,col.names = F,sep='\t')
write.table(data.frame(granges(joindeserts))[,1:3],'deserts_joined.bed',quote = F,row.names = F,col.names = F,sep='\t')
#----------------------------------------------------------------#
# Figure 1 B-C                                                   #
#----------------------------------------------------------------#
# ggplot theme
ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
theme_set(ggplot2_blanktheme)

# function for plotting signal across chromosomes
plot_signal <- function(sample,chrnum,tile) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = 'sacCer3')
  sample_ex <- sort(GenomeInfoDb::sortSeqlevels(sample))
  GenomeInfoDb::seqlengths(sample_ex) <- GenomeInfoDb::seqlengths(genome_info)
  bins_ex <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sample_ex),
                                       tilewidth=tile,
                                       cut.last.tile.in.chrom=TRUE)
  score_ex <- GenomicRanges::coverage(sample_ex, weight="score")
  bins_ex <- GenomicRanges::binnedAverage(bins_ex, score_ex, "binned_score")
  bins_ex <- GenomeInfoDb::keepSeqlevels(bins_ex, paste0("chr",chrnum),pruning.mode="coarse")
  positions_ex <- bins_ex@ranges@start + floor(bins_ex@ranges@width / 2)
  df_ex <- data.frame(seqnames=paste0('chr',chrnum),position=positions_ex, signal=bins_ex$binned_score)
}

library(patchwork)

#----------------------------------------------------------------#
# Figure 1 B
Red1_WT = hwglabr2::import_bedGraph("Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_WTd = gendiv(Red1_WT)

max_score=max(Red1_WTd$score)
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
p=list()
for(i in chrs) {
  print(i)
  p[[i]] = plot_signal(Red1_WTd,i,1000)
}

clusters = rtracklayer::import.bed('clusters_joined.bed')

plist=list()
for(i in chrs){
  a = data.frame()
  temp=data.frame()
  temp = p[[i]]
  if (exists("cen_mid")) {
    rm(cen_mid) }
  cen_mid = genome_info[seqnames(genome_info)==paste0("chr",i)]
  cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
  cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
  a <- ggplot(temp,aes(position,signal)) +
    geom_line(position='identity') +
    ylab(i) + ylim(-0.65,max_score) + xlim(0,1531933)
  a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                      size = 1.6, colour = 'green',shape=20)
  clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
  a <- a + geom_segment(clusterregion, size = 6,alpha = 0.6,
                        mapping=aes(x = start, y = max_score/2, xend = end, yend = max_score/2, colour = "segment"))
  plist[[i]] <- a
}

wrap_plots(plist,nrow = 16)

#----------------------------------------------------------------#
# Figure 1 C
Red1_rec8 <- hwglabr2::import_bedGraph("Red1-rec8D-39-62-193-90-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_rec8d = gendiv(Red1_rec8)

max_score=max(Red1_rec8d$score)
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
p=list()
for(i in chrs) {
  print(i)
  p[[i]] = plot_signal(Red1_rec8d,i,1000)
}

clusters = rtracklayer::import.bed('clusters_joined.bed')

plist=list()
for(i in chrs){
  a = data.frame()
  temp=data.frame()
  temp = p[[i]]
  if (exists("cen_mid")) {
    rm(cen_mid) }
  cen_mid = genome_info[seqnames(genome_info)==paste0("chr",i)]
  cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
  cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
  a <- ggplot(temp,aes(position,signal)) +
    geom_line(position='identity') +
    ylab(i) + ylim(-0.65,max_score) + xlim(0,1531933)
  a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                      size = 1.6, colour = 'green',shape=20)
  clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
  a <- a + geom_segment(clusterregion, size = 6,alpha = 0.6,
                        mapping=aes(x = start, y = max_score/2, xend = end, yend = max_score/2, colour = "segment"))
  plist[[i]] <- a
}

wrap_plots(plist,nrow = 16)

#----------------------------------------------------------------#
# Figure 1 D                                                     #
#----------------------------------------------------------------#
Red1_WT = hwglabr2::import_bedGraph("Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_WTd = gendiv(Red1_WT)

# Function to calculate average in regions
getmean <- function(regions, signalfile) {
  regionsignals = vector()
  for(i in 1:length(regions)){
    regionsignals[i] = regioneR::meanInRegions(A=regions[i], x=signalfile)
  }
  return(regionsignals)
}
# calcuate Red1 averages in clusters and deserts
Red1_WTd_in_cluster = getmean(clusters,Red1_WTd)
Red1_WTd_in_desert = getmean(deserts,Red1_WTd)

cluster_Red1_WTd_df = data.frame(region = 'cluster',signal = Red1_WTd_in_cluster)
desert_Red1_WTd_df = data.frame(region = 'desert',signal = Red1_WTd_in_desert)
Red1_all = rbind(cluster_Red1_WTd_df,desert_Red1_WTd_df)

# Plot Red1 averages in clusters and deserts
library(ggplot2)
p <- ggplot(Red1_all, aes(x=region, y=signal, fill = region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") +
  ylab("Red1 signal") + xlab("") + geom_boxplot(fill=NA)
p
t.test(Red1_WTd_in_cluster,Red1_WTd_in_desert) #p-value < 2.2e-16

#----------------------------------------------------------------#
# Figure 1 E                                                     #
#----------------------------------------------------------------#
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

flank <- floor(GenomicRanges::width(clusters)/2)
GenomicRanges::start(clusters) <- GenomicRanges::start(clusters) - flank
GenomicRanges::end(clusters) <- GenomicRanges::end(clusters) + flank
clusters <- GenomicRanges::trim(clusters)

flank <- floor(GenomicRanges::width(deserts)/2)
GenomicRanges::start(deserts) <- GenomicRanges::start(deserts) - flank
GenomicRanges::end(deserts) <- GenomicRanges::end(deserts) + flank
deserts <- GenomicRanges::trim(deserts)

Red1_WT = import_bedGraph("Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1_WTd = gendiv(Red1_WT)

# get signal in meta-regions
n_windows <- floor(1000 / 10)
Red1_clustermat <- EnrichedHeatmap::normalizeToMatrix(Red1_WTd, clusters,
                                                      value_column="score",
                                                      mean_mode="absolute",
                                                      extend=0, k=n_windows,
                                                      empty_value=NA, smooth=FALSE,
                                                      target_ratio=1)
Red1_clustermat_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Red1_clustermat,
                                                        ci=0.95, rep_bootstrap=1000,
                                                        na_rm=TRUE)
Red1_desertmat <- EnrichedHeatmap::normalizeToMatrix(Red1_WTd, deserts,
                                                     value_column="score",
                                                     mean_mode="absolute",
                                                     extend=0, k=n_windows,
                                                     empty_value=NA, smooth=FALSE,
                                                     target_ratio=1)
Red1_desertmat_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Red1_desertmat,
                                                       ci=0.95, rep_bootstrap=1000,
                                                       na_rm=TRUE)

group1s_gg <- data.frame(Data="cluster",Position=seq(1, n_windows), Red1_clustermat_metaORF)
group2s_gg <- data.frame(Data="desert",Position=seq(1, n_windows), Red1_desertmat_metaORF)
groups_metaregion = rbind(group1s_gg,group2s_gg)

# plot meta region signals
quarter1 <- round(n_windows / 4)
quarter2 <- quarter1 * 2
attr(Red1_clustermat, "upstream_index") <- 1:quarter1
attr(Red1_clustermat, "target_index") <- (quarter1 + 1):quarter2
attr(Red1_clustermat, "downstream_index") <- (quarter2 + 1):n_windows
attr(Red1_desertmat, "upstream_index") <- 1:quarter1
attr(Red1_desertmat, "target_index") <- (quarter1 + 1):quarter2
attr(Red1_desertmat, "downstream_index") <- (quarter2 + 1):n_windows

p <- ggplot(groups_metaregion, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(25,75), lty = 3) +
  scale_x_continuous(breaks = c(25,75),
                     labels = c('upstream','downstream'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p

#----------------------------------------------------------------#
# Figure 1 F                                                     #
#----------------------------------------------------------------#
# function to calculate peak density
peakdensity = function(peaksdata) {
  clusters = rtracklayer::import.bed('clusters_joined.bed')
  deserts = rtracklayer::import.bed('deserts_joined.bed')

  cluster_peaks = vector()
  for(i in 1:length(clusters)) {
    hits <- findOverlaps(query = clusters[i],subject = peaksdata)
    overlaps <- pintersect(clusters[i][queryHits(hits)], peaksdata[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(peaksdata[subjectHits(hits)])
    hits <- hits[percentOverlap > 0.5]
    cluster_peaks[i] <- length(peaksdata[subjectHits(hits)])/width(clusters[i])*1000
    rm(hits);rm(overlaps)
  }
  desert_peaks = vector()
  for(i in 1:length(deserts)) {
    hits <- findOverlaps(query = deserts[i],subject = peaksdata)
    overlaps <- pintersect(deserts[i][queryHits(hits)], peaksdata[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(peaksdata[subjectHits(hits)])
    hits <- hits[percentOverlap > 0.5]
    desert_peaks[i] <- length(peaksdata[subjectHits(hits)])/width(deserts[i])*1000
    rm(hits);rm(overlaps)
  }
  clusterpeaks = data.frame(region="cluster",density = cluster_peaks)
  desertpeaks = data.frame(region="desert",density = desert_peaks)
  allpeaks = rbind(clusterpeaks,desertpeaks)
  return(allpeaks)
}

# Calculate density for Red1
Red1peaks = read.table('Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
Red1peaks = GRanges(Red1peaks$V1,IRanges(Red1peaks$V2,Red1peaks$V3))

allpeaksRed1 = peakdensity(Red1peaks)
t.test(allpeaksRed1[which(allpeaksRed1$region=='cluster'),'density'],allpeaksRed1[which(allpeaksRed1$region=='desert'),'density']) #p-value = 7.371e-07
wilcox.test(allpeaksRed1[which(allpeaksRed1$region=='cluster'),'density'],allpeaksRed1[which(allpeaksRed1$region=='desert'),'density']) #p-value = 1.945e-07

# Calculate density for Hop1
Hop1peaks = read.table('Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
Hop1peaks = GRanges(Hop1peaks$V1,IRanges(Hop1peaks$V2,Hop1peaks$V3))

allpeaksHop1 = peakdensity(Hop1peaks)
t.test(allpeaksHop1[which(allpeaksHop1$region=='cluster'),'density'],allpeaksHop1[which(allpeaksHop1$region=='desert'),'density']) #p-value = 5.706e-06
wilcox.test(allpeaksHop1[which(allpeaksHop1$region=='cluster'),'density'],allpeaksHop1[which(allpeaksHop1$region=='desert'),'density']) #p-value = 2.437e-06

# Calculate density for Hop1
Rec8peaks = read.table('Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
Rec8peaks = GRanges(Rec8peaks$V1,IRanges(Rec8peaks$V2,Rec8peaks$V3))

allpeaksRec8 = peakdensity(Rec8peaks)
t.test(allpeaksRec8[which(allpeaksRec8$region=='cluster'),'density'],allpeaksRec8[which(allpeaksRec8$region=='desert'),'density']) #p-value = 0.0001465
wilcox.test(allpeaksRec8[which(allpeaksRec8$region=='cluster'),'density'],allpeaksRec8[which(allpeaksRec8$region=='desert'),'density']) #p-value = 9.813e-06

# plot all densities
combpeaksRed1 = data.frame(Region = paste0('Red1',allpeaksRed1$region),density = allpeaksRed1$density)
combpeaksHop1 = data.frame(Region = paste0('Hop1',allpeaksHop1$region),density = allpeaksHop1$density)
combpeaksRec8 = data.frame(Region = paste0('Rec8',allpeaksRec8$region),density = allpeaksRec8$density)
combpeaks = rbind(combpeaksRed1,combpeaksHop1,combpeaksRec8)
ggplot(combpeaks, aes(Region, density, group=Region, fill = Region)) +
  geom_violin()+ geom_boxplot(fill=NA)

#----------------------------------------------------------------#
# Figure 1 G                                                     #
#----------------------------------------------------------------#
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
# G: Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# H: Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# I: Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# J: AH8867B-379-767-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# K: AH6408I-144-183-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# L: AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz
# M: Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz

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
