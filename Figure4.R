#----------------------------------------------------------------#
# Figure 4A                                                      #
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
# Figure 4B                                                      #
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
# Figure 4 C-D                                                   #
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
library(caret)
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
# Figure 4C                                                      #
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
p <- ggplot(tempm, aes(x=codingdensity, y=clusterdensity)) + geom_dotplot(binwidth = 1, y=2403) +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

q <- ggplot(tempm, aes(x=codingdensity, y=clusterdensity)) + geom_dotplot(binwidth = 1, y=2403) 

library(beeswarm)
clusterTemp <- subset(tempm, tempm$clusterdensity > 0)
desertTemp <- subset(tempm, tempm$clusterdensity < 0.5)
clusterTemp$region = "Island"
desertTemp$region = "Desert"
Toget <- rbind(clusterTemp, desertTemp)


beeswarm(codingdensity ~ region, data = Toget, pch = 16, 
         col = c("deepskyblue1", "red"), xlab = "Coding Density", ylab = "Region",
         corral = "wrap", vertical = FALSE)

b <- ggplot(Toget, aes(x=region, y=codingdensity)) + geom_dotplot(binwidth = 1, y=)

r <- ggplot(desertTemp, aes(codingdensity)) + ylim(0,75)
r <- r+ geom_histogram(bins = 100, color = "black", fill = "deepskyblue1")
r <- r+ geom_vline(aes(xintercept=mean(codingdensity)), linetype="dashed")
r <- r + ggtitle("Desert Regions")

k <- ggplot(clusterTemp, aes(codingdensity))+ ylim(0,75)+ geom_histogram(bins = 100, color = "black", fill = "red")

k <- k+ geom_vline(aes(xintercept=mean(codingdensity)), linetype="dashed")
k <- k + ggtitle("Island Regions")

library(patchwork)

Plot <- k / r