#---------------------------------------------------------------#
# Axis clusters                                                 #
# Supplemental Figure 3 code                                    #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)

#---------------------------------------------------------------#
# Fig S3 A                                                      #
#---------------------------------------------------------------#
peakdensity = function(peaksdata) {
  clusters = rtracklayer::import.bed('clusters_joined.bed')
  deserts = rtracklayer::import.bed('deserts_joined.bed')

  cluster_peaks = vector()
  for(i in 1:length(clusters)) {
    #print(clusters[i])
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


hotspots = get_dsb_hotspots('SK1Yue')
allspo11oligos = peakdensity(hotspots)

ggplot(allspo11oligos, aes(region, density, group=region, fill = region)) +
  geom_violin()+ geom_boxplot(fill=NA)
hist(allspo11oligos$density)
wilcox.test(allspo11oligos[which(allspo11oligos$region=='cluster'),'density'],allspo11oligos[which(allspo11oligos$region=='desert'),'density']) #p-value = 0.389

#---------------------------------------------------------------#
# Fig S3 B                                                      #
#---------------------------------------------------------------#
spo11oligos = rtracklayer::import.bedGraph('Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg')
genAvg <- hwglabr2::average_chr_signal(spo11oligos)$genome_avrg
spo11oligos$score <- spo11oligos$score/genAvg

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

matdesertspo11oligos <- EnrichedHeatmap::normalizeToMatrix(spo11oligos, hotspots_all[which(hotspots_all$class=='desert')], value_column = "score",
                                                  extend = 300, mean_mode = "weighted", w = 1,empty_value=0)
matclusterspo11oligos <- EnrichedHeatmap::normalizeToMatrix(spo11oligos, hotspots_all[which(hotspots_all$class=='cluster')], value_column = "score",
                                                   extend = 300, mean_mode = "weighted", w = 1,empty_value=0)
###
auc1data <- c(mean(rowSums(data.frame(matdesertspo11oligos)),na.rm=T),
              mean(rowSums(data.frame(matclusterspo11oligos)),na.rm=T))
auc1sd <- c(sd(rowSums(data.frame(matdesertspo11oligos)),na.rm=T)/sqrt(length(rowSums(data.frame(matdesertspo11oligos)))),
            sd(rowSums(data.frame(matclusterspo11oligos)),na.rm=T)/sqrt(length(rowSums(data.frame(matclusterspo11oligos)))))
cvdgroups <- c('spo11oligos_desert','spo11oligos_cluster')
auc1df <- data.frame(cvdgroups,auc1data,auc1sd)
auc1df$cvdgroups <- factor(auc1df$cvdgroups, levels=unique(auc1df$cvdgroups))
auc1df$cvdgroups <- factor(auc1df$cvdgroups, levels=c('spo11oligos_desert','spo11oligos_cluster'))

ggplot(auc1df,aes(x=cvdgroups,y=auc1data,fill=cvdgroups,width=0.8)) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_errorbar(aes(ymin=auc1data-auc1sd, ymax=auc1data+auc1sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  labs(title = '', x = '', y = 'Relative Spo11 oligo signal')

wilcox.test(rowSums(data.frame(matdesertspo11oligos)),rowSums(data.frame(matclusterspo11oligos)))  #p-value = 0.04847

#---------------------------------------------------------------#
# Fig S3 C                                                      #
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

bedgraphs <- 'H3-WT-44-105-445-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz'

dir.create("hotspots_pdf")
# function to plot average signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][4],split="_")[[1]][1]

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
# Fig S3 D                                                      #
#---------------------------------------------------------------#
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

gff <- get_gff("SK1Yue")
gff <- gff[gff$type=='gene']
transcription <- read.csv('2017.06.16_SK1Yue_EdgeR_tpm.csv')
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]
gff_txn <- GRanges(seqnames = gff_txn$seqnames,IRanges(gff_txn$start,gff_txn$end),
                   score=gff_txn$AH119_3h,ID=gff_txn$ID)

hits <- findOverlaps(clusters,gff_txn)
overlaps <- pintersect(clusters[queryHits(hits)], gff_txn[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff_txn[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
cluster_gff_txn <- gff_txn[subjectHits(hits)]
rm(hits);rm(overlaps)

hits <- findOverlaps(deserts,gff_txn)
overlaps <- pintersect(deserts[queryHits(hits)], gff_txn[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff_txn[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
deserts_gff_txn <- gff_txn[subjectHits(hits)]

ctxn <- data.frame(Type="cluster",txn=cluster_gff_txn$score)
dtxn <- data.frame(Type="desert",txn=deserts_gff_txn$score)
alltxn <- rbind(ctxn,dtxn)
p <- ggplot(alltxn, aes(x=Type,y=txn, fill = Type)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") + ylim(0,400) +
  ylab("Txn") + xlab("") + geom_boxplot(fill=NA)
p
wilcox.test(deserts_gff_txn$score,cluster_gff_txn$score) #p-value = 3.042e-06
