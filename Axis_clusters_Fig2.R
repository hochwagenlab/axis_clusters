#---------------------------------------------------------------#
# Axis clusters                                                 #
# Figure 2 code                                                 #
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
# Figure 2 A                                                     #
#----------------------------------------------------------------#

clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
gff <- hwglabr2::get_gff("SK1Yue")
gff <- gff[gff$type=='gene']

# cluster genes
hits <- findOverlaps(clusters,gff)
overlaps <- pintersect(clusters[queryHits(hits)], gff[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
cluster_gff <- gff[subjectHits(hits)]
rm(hits);rm(overlaps)

# desert genes
hits <- findOverlaps(deserts,gff)
overlaps <- pintersect(deserts[queryHits(hits)], gff[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(gff[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
deserts_gff <- gff[subjectHits(hits)]
rm(hits);rm(overlaps)

cluster_genesize <- data.frame(Type = 'gene_cluster',Region="cluster",width=width(cluster_gff))
desert_genesize <- data.frame(Type = 'gene_desert',Region="desert",width=width(deserts_gff))
allsize <- rbind(cluster_genesize,desert_genesize)
p <- ggplot(allsize, aes(x=Region,y=width, group=Region, fill = Region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") +
  ylab("Size (bp)") + xlab("") + geom_boxplot(fill=NA)
p
wilcox.test(width(cluster_gff),width(deserts_gff)) #p-value = 4.293e-13

#----------------------------------------------------------------#
# Figure 2 B                                                     #
#----------------------------------------------------------------#
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

intergen <- get_intergenic_regions("SK1Yue",as_gr = T)

# cluster intergens
hits <- findOverlaps(clusters,intergen)
overlaps <- pintersect(clusters[queryHits(hits)], intergen[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(intergen[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
cluster_intergen <- intergen[subjectHits(hits)]
rm(hits);rm(overlaps)

# desert intergens
hits <- findOverlaps(deserts,intergen)
overlaps <- pintersect(deserts[queryHits(hits)], intergen[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(intergen[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
deserts_intergen <- intergen[subjectHits(hits)]

# divide into type of intergen
div_desert = deserts_intergen[deserts_intergen$type == 'divergent']
tan_desert = deserts_intergen[deserts_intergen$type == 'tandem']
conv_desert = deserts_intergen[deserts_intergen$type == 'convergent']
div_cluster = cluster_intergen[cluster_intergen$type == 'divergent']
tan_cluster = cluster_intergen[cluster_intergen$type == 'tandem']
conv_cluster = cluster_intergen[cluster_intergen$type == 'convergent']

div_desert_df = data.frame(intergen='divergent-desert',region = 'desert',width = width(div_desert))
conv_desert_df = data.frame(intergen='convergent-desert',region = 'desert',width = width(conv_desert))
tan_desert_df = data.frame(intergen='tandem-desert',region = 'desert',width = width(tan_desert))
div_cluster_df = data.frame(intergen='divergent-cluster',region = 'cluster',width = width(div_cluster))
conv_cluster_df = data.frame(intergen='convergent-cluster',region = 'cluster',width = width(conv_cluster))
tan_cluster_df = data.frame(intergen='tandem-cluster',region = 'cluster',width = width(tan_cluster))
intergen_all = rbind(conv_cluster_df,conv_desert_df,tan_cluster_df,tan_desert_df,div_cluster_df,div_desert_df)

p <- ggplot(intergen_all, aes(x=intergen, y=width, fill = region)) +
  geom_violin()+scale_fill_brewer(palette="Dark2") + ylim(0,2000) +
  ylab("Width of intergens") + xlab("") + geom_boxplot(fill=NA)
p

#----------------------------------------------------------------#
# Figure 2 C-D                                                   #
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
tempdf = data.frame(temp)
tempm = tempdf[,c('seqnames','clusterdensity','codingdensity')]

# Plot feature
featurePlot(x = tempm[, 'codingdensity'],
            y = tempm$clusterdensity,
            plot = "scatter",
            type = c("p", "smooth"))
ggplot(tempm, aes(x=codingdensity, y=clusterdensity)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

# logistic regression test
library(dplyr)
temp$clusterdensity = round(temp$clusterdensity)

default_idx = createDataPartition(tempm$clusterdensity, p = 0.8, list = FALSE)
default_trn = tempm[default_idx, ]
default_tst = tempm[-default_idx, ]

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
#Area under the curve: 0.6919

#----------------------------------------------------------------#
# Figure 2 E                                                     #
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


# Concert window_size to n of windows
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
