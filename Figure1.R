#---------------------------------------------------------------#
# Axis Islands -- Islands = clusters in code                    #
# Figure 1 code                                                 #
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
# Figure 1A (for all chromosomes)                                #
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
# Figure 1B                                                      #
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
# Figure 1C                                                      #
#----------------------------------------------------------------#
Hop1_WT = hwglabr2::import_bedGraph("Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Hop1_WTd = gendiv(Hop1_WT)

# Function to calculate average in regions
getmean <- function(regions, signalfile) {
  regionsignals = vector()
  for(i in 1:length(regions)){
    regionsignals[i] = regioneR::meanInRegions(A=regions[i], x=signalfile)
  }
  return(regionsignals)
}
# calcuate Hop1 averages in clusters and deserts
Hop1_WTd_in_cluster = getmean(clusters,Hop1_WTd)
Hop1_WTd_in_desert = getmean(deserts,Hop1_WTd)

cluster_Hop1_WTd_df = data.frame(region = 'cluster',signal = Hop1_WTd_in_cluster)
desert_Hop1_WTd_df = data.frame(region = 'desert',signal = Hop1_WTd_in_desert)
Hop1_all = rbind(cluster_Hop1_WTd_df,desert_Hop1_WTd_df)

# Plot Hop1 averages in clusters and deserts
library(ggplot2)
p <- ggplot(Hop1_all, aes(x=region, y=signal, fill = region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") +
  ylab("Hop1 signal") + xlab("") + geom_boxplot(fill=NA)
p
t.test(Hop1_WTd_in_cluster,Hop1_WTd_in_desert) #p-value < 2.2e-16


#----------------------------------------------------------------#
# Figure 1D E                                                    #
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
# D: Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
#    Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# E: Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
#    AH8867B-379-767-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz


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

#Figure 1F


ggplot2_theme <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=10),
        plot.subtitle=element_text(hjust=0.5, size=10),
        axis.text=element_text(colour='black'),
        axis.ticks=element_line(colour='black'))
theme_set(ggplot2_theme)

clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
gff <- get_gff("SK1Yue")
gff <- gff[gff$type=='gene']
Red1 <- import_bedGraph('Red1-wildtype-71-34-199-29-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Red1 = gendiv(Red1)

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

convClus <- read.csv('dfConv_cluster.csv')
convClus <- toGRanges(convClus)
convCls_gff <- subsetByOverlaps(cluster_gff, convClus)

convDest <- read.csv('/Users/darmokandjalad/Documents/HI-Scripts_Analysis/IslandPaper/dfConv_desert.csv')
convDest <- toGRanges(convDest)
convDes_gff <- subsetByOverlaps(deserts_gff, convDest)


sig_at_convClut <- hwglabr2::signal_at_orf2(signal_data=Red1, gff=convCls_gff,
                                            write_to_file=FALSE)
sig_at_convCluts <- hwglabr2::signal_mean_and_ci(signal_data=sig_at_convClut,
                                                 ci=0.95, rep_bootstrap=1000,
                                                 na_rm=TRUE)
sig_at_convDst <- hwglabr2::signal_at_orf2(signal_data=Red1, gff=convDes_gff,
                                           write_to_file=FALSE)
sig_at_convDsts <- hwglabr2::signal_mean_and_ci(signal_data=sig_at_convDst,
                                                ci=0.95, rep_bootstrap=1000,
                                                na_rm=TRUE)

Clust_gg <- data.frame(Data="cIslands",Position=seq(1, 1000), sig_at_convCluts)
Dests_gg <- data.frame(Data="Deserts",Position=seq(1, 1000), sig_at_convDsts)

allgroups <- rbind(Clust_gg,Dests_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p

p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()


###now try all of them instead of just convergent
sig_at_convClut <- hwglabr2::signal_at_orf2(signal_data=Red1, gff=cluster_gff,
                                            write_to_file=FALSE)
sig_at_convCluts <- hwglabr2::signal_mean_and_ci(signal_data=sig_at_convClut,
                                                 ci=0.95, rep_bootstrap=1000,
                                                 na_rm=TRUE)
sig_at_convDst <- hwglabr2::signal_at_orf2(signal_data=Red1, gff=deserts_gff,
                                           write_to_file=FALSE)
sig_at_convDsts <- hwglabr2::signal_mean_and_ci(signal_data=sig_at_convDst,
                                                ci=0.95, rep_bootstrap=1000,
                                                na_rm=TRUE)

Clust_gg <- data.frame(Data="Islands",Position=seq(1, 1000), sig_at_convCluts)
Dests_gg <- data.frame(Data="Deserts",Position=seq(1, 1000), sig_at_convDsts)

allgroups <- rbind(Clust_gg,Dests_gg)

# Set up the plot
b <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
b <- b + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
b

p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()




