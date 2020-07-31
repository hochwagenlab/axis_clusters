#---------------------------------------------------------------#
# Axis clusters                                                 #
# Supplemental Figure 2 code                                    #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)

#----------------------------------------------------------------#
# Figure S2 A-B                                                  #
#----------------------------------------------------------------#
# Make 5kb tiles, drop the last tile per chromosomes (because not = 5000)
tilesize=5000
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_tiles <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(genome_info),
                                          tilewidth=tilesize,
                                          cut.last.tile.in.chrom=TRUE)
genome_tilemidpoints = genome_tiles[width(genome_tiles)==tilesize]
midpoint <- floor(width(genome_tilemidpoints) / 2)
start(genome_tilemidpoints) <- start(genome_tilemidpoints) + midpoint
end(genome_tilemidpoints) <- start(genome_tilemidpoints)

# cluster density
clusters = rtracklayer::import.bed('clusters_joined.bed')
score(clusters) = 1.0
clusterdensity <- normalizeToMatrix(signal = clusters, target = genome_tilemidpoints, value_column = "score",
                                    extend = c(tilesize/2-1,tilesize/2), mean_mode = "coverage", w = 1,background = NA)
mcols(genome_tilemidpoints)['clusterdensity'] <- rowSums(data.frame(clusterdensity),na.rm = T)/ncol(data.frame(clusterdensity))

# coding density
gff = hwglabr2::get_gff('SK1Yue')
gff = granges(gff[gff$type=='gene'])
score(gff) = 1
coding <- normalizeToMatrix(gff, genome_tilemidpoints, value_column = "score",
                            extend = c(tilesize/2-1,tilesize/2), mean_mode = "coverage",
                            w = 1,background=NA)
mcols(genome_tilemidpoints)['count'] <- rowSums(data.frame(coding),na.rm = T)
mcols(genome_tilemidpoints)['codingdensity'] <- rowSums(data.frame(coding),na.rm = T)/ncol(data.frame(coding))

temp = genome_tilemidpoints
temp$clusterdensity = round(temp$clusterdensity)

# remove small chromosomes and test chr
testchr = 'chrIX'
largechrclusters = temp[seqnames(temp)!='chrI']
largechrclusters = largechrclusters[seqnames(largechrclusters)!='chrIII']
largechrclusters = largechrclusters[seqnames(largechrclusters)!='chrVI']
largechrclustersdf = data.frame(largechrclusters)[,c('clusterdensity','codingdensity')]


# Plot feature
ggplot(largechrclustersdf, aes(x=codingdensity, y=clusterdensity)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

# Logistic regression model
default_idx = createDataPartition(largechrclustersdf$clusterdensity, p = 0.8, list = FALSE)
default_trn = largechrclustersdf[default_idx, ]
default_tst = largechrclustersdf[-default_idx, ]

glm_mod = train(
  as.factor(clusterdensity) ~ codingdensity,
  data = default_trn,
  method = "glm",
  family = "binomial"
)

pred.glmModel <- predict(glm_mod, newdata=default_tst, type="prob")
library(pROC)
roc.glmModel <- roc(default_tst$clusterdensity, pred.glmModel[,'1'])
plot(roc.glmModel,col='orange')
auc.glmModel <- auc(roc.glmModel)
# Area under the curve: 0.7252

tempdf = data.frame(temp)
tempm = tempdf[,c('seqnames','clusterdensity','codingdensity')]

#---------------------------------------------------------------#
# Fig S2 C-E, G-H                                               #
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
# C: AH7797_MNaseSeq_0h-SK1Yue-PM_B3-MACS2_treat_pileup.bdg
# D: WT_Hu2015_MNaseSeq-SK1Yue-PM_B3-MACS2_treat_pileup.bdg
# E: H4K44R_Hu2015_MNaseSeq-SK1Yue-PM_B3-MACS2_treat_pileup.bdg
# G: Nucleosome_reps-SK1-MACS2_treat_pileup.bdg
# H: AH9003_MNaseSeq_3h-SK1Yue-PM_B3-MACS2_treat_pileup.bdg

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
# Fig S2 F                                                      #
#---------------------------------------------------------------#

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
bedgraphs <- "AH9003_MNaseSeq_3h-SK1Yue-PM_B3-MACS2_treat_pileup.bdg"

# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]

  Red1_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Red1_bg)$genome_avrg
  Red1_bg$score <- Red1_bg$score/genAvg

  mat1 <- normalizeToMatrix(Red1_bg, hotspots_all, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

  col_fun <- circlize::colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))

  pdf(paste0("hotspots_pdf/MNase_rtf1_around_hotspots_heatmap_sort.pdf"))
  print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c('green','red')),
                                                                                 show_error = TRUE,pos_line=FALSE)),
                        row_title_rot = 0,
                        axis_name = c("-1 kb", "hotspots", "1 kb"),
                        row_order = 1:length(hotspots_all),
                        split=hotspots_all$class,
                        column_title ="MNase in rtf1"))
  dev.off()
}

#---------------------------------------------------------------#
# Fig S2I-K                                                     #
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
# F: AH9003-Red1chip-191219-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# G: AH8584-Red1-473-574-784-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz
# H: AH8583-Red1-475-576-781-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz

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
