library(Seurat)
library(CountClust)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)

# CITATION: Most code below was adapted from or inspired by 
# the Broad single-cell workshop tutorial (2020): 
# https://broadinstitute.github.io/2020_scWorkshop/feature-selection-and-cluster-analysis.html#probabilistic-lda-clustering

## For reference
# topic12/X12: IFN response
# topic14/X14: Secretory
# topic15/X15: APC

#######################################################
########## Fit topics in Foxn1 microglia ##############
#######################################################

# read in full Foxn1 myeloid Seurat object
myeloid.unintegrated <- readRDS("~/Desktop/Microglia_Project_Current/myeloid.unintegrated.cr3.rds")

# subset to microglia clusters
microglia_subset <- subset(myeloid.unintegrated, seurat_clusters %in% c("0","1","2"))

# ensure assay is set to RNA
DefaultAssay(microglia_subset) <- "RNA"

# extract counts matrix
microglia_counts <- as.matrix(GetAssayData(microglia_subset, slot = "counts", assay = "RNA"))

# fit 15 topics with a tolerance of 10
microglia_pbmc_fitGoM_K15 <- FitGoM(t(microglia_counts), K=15, tol=10)

# extract topic weights & add to metadata
omega <- data.frame(microglia_pbmc_fitGoM_K15$fit$omega)
microglia_subset  <- AddMetaData(microglia_subset ,omega)

# extract gene weights
theta_mat <- microglia_pbmc_fitGoM_K15$fit$theta

################################################
################ Figure scripts ################
################################################

###############################################
####### Figure 2 - Foxn1 topic summary ########
###############################################

########## top topic gene marker heatmap ##########

# set up empty lists & vectors
topic_df <- list()
gene_list <- c()

# extract top 15 gene markers for each topic &
# save list of top genes & topic data frame list
for(i in c(12,14,15)){
  top_features <- ExtractTopFeatures(theta_mat,
                                     top_features = 15,
                                     method = "poisson",
                                     options="min",
                                     shared = TRUE)
  
  topic_df[[i]] <- data.frame(scores = top_features$scores[i,],
                              genes = rownames(theta_mat)[top_features$indices[i,]],
                              topic = rep(paste0("topic",i),length(top_features$scores[i,])))
  
  gene_list <- c(gene_list, rownames(theta_mat)[top_features$indices[i,]])
}

# plot marker heatmap
pheatmap(
  # extract top genes from original topic/gene weights matrix
  t(microglia_pbmc_fitGoM_K15$fit$theta[unique(gene_list),c(12,14,15)]), 
  # scale gene weights across the 3 topics in the plot
  scale = "column", 
  # cluster topics & genes
  cluster_rows = T, cluster_cols = T, 
  # set color palette
  color = colorRampPalette(c("white","white",brewer.pal(n = 5, name ="Greys")))(100))


######### 3 topic overlay t-SNE plot ###########

# add tsne embedding to @metadata slot
microglia_subset <- AddMetaData(microglia_subset, 
                                data.frame(microglia_subset@reductions$tsne@cell.embeddings))

# extract & melt topic/tsne info for plot
plot_df <- microglia_subset@meta.data[,c("tSNE_1","tSNE_2","X12","X14","X15")] %>% 
  gather(X12:X15, key="topic",value = "score")

### full plot statement
ggplot() + 
  # plot all cells in transparent gray
  geom_point(data=plot_df, aes(x=tSNE_1, y=tSNE_2), color="gray90",alpha=0.1, size=0.5) + 
  # add another layer of points that's subset to cells that have a topic weight 
  # greater than .1 in at least one of the three topics of interest
  # and enforce that the transparency will be inversely related 
  # to how high the topic score is (alpha = score)
  geom_point(data=subset(plot_df,score > 0.1), 
             aes(x=tSNE_1, y=tSNE_2, color = topic, alpha=score), size=0.5) + 
  # set custom colors for each topic
  scale_color_manual(values = c("X14" = "lightpink2", 
                                "X12" = "#1E87E5FF", 
                                "X15" = "yellow3")) + 
  # back background white & remove t-SNE axes
  theme_classic() + NoAxes()

########## individual topic FeaturePlots ##########

FeaturePlot(microglia_subset,"X12", order=TRUE, pt.size = 0.1, min.cutoff = 0) & 
  NoAxes() & ggtitle("") & scale_color_gradientn(colors = c(pal_material("grey")(10)[3:10]))
FeaturePlot(microglia_subset,"X14", order=TRUE, pt.size = 0.1, min.cutoff = 0) & 
  NoAxes() & ggtitle("")  & scale_color_gradientn(colors = c(pal_material("grey")(10)[3:10]))
FeaturePlot(microglia_subset,"X15", order=TRUE, pt.size = 0.1, min.cutoff = 0) & 
  NoAxes() & ggtitle("")  & scale_color_gradientn(colors = c(pal_material("grey")(10)[3:10]))

###################################################
####### Figure S2 - enriched topic barplot ########
###################################################

# adapt summary data of relative topic enrichments & 
# logic for truly 'Met-enriched' sets
plot_df <- 
  unique(microglia_subset@meta.data[,c("orig.ident", paste0("X",1:15))] %>% 
           # set topics to be a row
           gather(X1:X15, key="topic",value = "assignment") %>% 
           # add back in status as a variable
           mutate(status = ifelse(orig.ident %in% c("C1","C2","C3"), 
                                  "Control","Metastatic")) %>% 
           group_by(status, topic) %>% 
           mutate(avg_topic_by_status = mean(assignment)) %>% 
           group_by(orig.ident, topic) %>% 
           mutate(avg_topic_by_mouse = mean(assignment)) %>% 
           # expand out status and the average value for each mouse
           spread(status, avg_topic_by_mouse) %>% 
           group_by(topic) %>% 
           # determine if a topic is "met_enriched" based on the
           # highest ctrl value and lowest met value across mice
           mutate(met_enriched = 
                    ifelse(max(Control, na.rm = TRUE) < min(Metastatic, na.rm=TRUE),
                           "Met_enriched","Non_enriched")) %>% 
           # move back to melted data
           gather(Control:Metastatic, key = "status", value = "avg_topic_by_mouse") %>%
           # remove NAs added by gather
           filter(!is.na(avg_topic_by_mouse)) %>%
           dplyr::select(-assignment)) %>%
  group_by(topic) %>%
  mutate(ctrl_topic_values = ifelse(status == "Control", avg_topic_by_mouse,NA)) %>%
  mutate(control_mean = mean(ctrl_topic_values, na.rm = TRUE)) %>%
  mutate(RQ = avg_topic_by_mouse - control_mean) 
# reorder topics to go from 1-15
plot_df$topic <- factor(plot_df$topic, levels = paste0("X",1:15))

# generate topic RQ barplot (relative to control mouse mean)
ggplot(plot_df, 
       aes(x = topic, y = RQ, 
           fill = orig.ident, 
           alpha = met_enriched)) + theme_bw() + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = c("C1" = "cadetblue1", 
                               "C2" = "cadetblue3", 
                               "C3" = "cadetblue4", 
                               "Met1" = "firebrick1", 
                               "Met2" = "firebrick3", 
                               "Met3" = "firebrick4")) + 
  ylab("Relative topic assignment") + 
  # make non-enriched topics transparent
  scale_alpha_manual(values = c("Non_enriched" = I(0.25), 
                                "Met_enriched" = I(1)))


############################################################
####### Figure S8 - Foxn1 topics in MITRG microglia ########
############################################################

# read in MITRG Seurat object
MITRG_myeloid.unintegrated <- readRDS("~/Desktop/human.mg.all.barcodes.rds")
# subset to microglia
MITRG_microglia_subset <- subset(MITRG_myeloid.unintegrated, 
                                 annotated_cluster %in% c("Homeostatic.MG",
                                                          "IFN.Responsive.MG"))


# extract top 25 marker genes for each Foxn1 topic
topic_df <- list()
gene_list <- c()
for(i in c(3,12,14,15)){
  top_features <- ExtractTopFeatures(theta_mat,
                                     top_features = 25,
                                     method = "poisson",
                                     options="min",
                                     shared = FALSE)
  
  topic_df[[i]] <- data.frame(scores = top_features$scores[i,],
                              genes = rownames(theta_mat)[top_features$indices[i,]])
  
  gene_list <- c(gene_list, rownames(theta_mat)[top_features$indices[i,]])
}


# gene symbol conversion function from mouse-to-human 
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x , 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2)
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

##### take top 25 unique, human-equivalent topic markers &
##### score MITRG microglia for them
topic3_genes_human <- convertMouseGeneList(as.character(topic_df[[3]]$genes))
topic3_genes_human <- unique(as.character(topic3_genes_human$HGNC.symbol))

MITRG_microglia_subset <- AddModuleScore(MITRG_microglia_subset, 
                                         features = list(topic3_genes_human),
                                         name = "Foxn1_Ribo_topic")

topic12_genes_human <- convertMouseGeneList(as.character(topic_df[[12]]$genes))
topic12_genes_human <- unique(as.character(topic12_genes_human$HGNC.symbol))
MITRG_microglia_subset <- AddModuleScore(MITRG_microglia_subset, 
                                         features = list(topic12_genes_human),
                                         name = "Foxn1_IFN_response_topic")

topic14_genes_human <- convertMouseGeneList(as.character(topic_df[[14]]$genes))
topic14_genes_human <- unique(as.character(topic14_genes_human$HGNC.symbol))
MITRG_microglia_subset <- AddModuleScore(MITRG_microglia_subset, 
                                         features = list(topic14_genes_human),
                                         name = "Foxn1_Secretory_topic")

topic15_genes_human <- convertMouseGeneList(as.character(topic_df[[15]]$genes))
topic15_genes_human <- unique(as.character(topic15_genes_human$HGNC.symbol))
MITRG_microglia_subset <- AddModuleScore(MITRG_microglia_subset, 
                                         features = list(topic15_genes_human),
                                         name = "Foxn1_APC_topic")

############### set up 3 topic score overlay plot #################
# (but instead of fit topic weights, now we use simple gene scores)

MITRG_microglia_subset <- AddMetaData(
  MITRG_microglia_subset, 
  data.frame(MITRG_microglia_subset@reductions$tsne@cell.embeddings))
plot_df <- MITRG_microglia_subset@meta.data[,c("tSNE_1","tSNE_2",
                                               "Foxn1_IFN_response_topic1",
                                               "Foxn1_Secretory_topic1",
                                               "Foxn1_APC_topic1")] %>% gather(Foxn1_IFN_response_topic1:Foxn1_APC_topic1, key="topic",value = "score") %>%
  mutate(zero_one_score = ifelse(score < 0, 0, score/max(score)))
#group_by(topic) %>%
#mutate(zero_one_score = (score - min(score))/(max(score) - min(score)))

ggplot() + geom_point(data=plot_df, aes(x=tSNE_1, y=tSNE_2), 
                      color="gray90",alpha=0.1, size=0.5) + 
  geom_point(data=subset(plot_df,zero_one_score > 0.25), 
             aes(x=tSNE_1, y=tSNE_2, color = topic, alpha=zero_one_score), size=0.5) + 
  scale_color_manual(values = c("Foxn1_Ribo_topic1" = "gray50", 
                                "Foxn1_Secretory_topic1" = "lightpink2", 
                                "Foxn1_IFN_response_topic1" = "#1E87E5FF", 
                                "Foxn1_APC_topic1" = "yellow3")) + 
  theme_classic() + NoAxes()

############### set up 3 topic score barplot summary #################

MITRG_plot_df <- 
  unique(MITRG_microglia_subset@meta.data[,c("multiseq.barcodes.adjusted",
                                             "Foxn1_IFN_response_topic1", 
                                             "Foxn1_Secretory_topic1", 
                                             "Foxn1_APC_topic1", 
                                             "Foxn1_Ribo_topic1")] %>% 
           # set topics to be a row
           gather(Foxn1_IFN_response_topic1:Foxn1_Ribo_topic1, key="topic",
                  value = "assignment") %>% 
           # add back in status as a variable
           mutate(status = ifelse(multiseq.barcodes.adjusted %in% c("Bar1","Bar2","Bar3"), 
                                  "Control","Metastatic")) %>% 
           group_by(status, topic) %>% 
           mutate(avg_topic_by_status = mean(assignment)) %>% 
           group_by(multiseq.barcodes.adjusted, topic) %>% 
           mutate(avg_topic_by_mouse = mean(assignment)) %>% 
           # expand out status and the average value for each mouse
           spread(status, avg_topic_by_mouse) %>% 
           group_by(topic) %>% 
           # determine if a topic is "met_enriched" based on the
           # highest ctrl value and lowest met value across mice
           mutate(met_enriched = 
                    ifelse(max(Control, na.rm = TRUE) < min(Metastatic, na.rm=TRUE),
                           "Met_enriched","Non_enriched")) %>% 
           # move back to melted data
           gather(Control:Metastatic, key = "status", value = "avg_topic_by_mouse") %>%
           # remove NAs added by gather
           filter(!is.na(avg_topic_by_mouse)) %>%
           dplyr::select(-assignment)) %>%
  group_by(topic) %>%
  mutate(ctrl_topic_values = ifelse(status == "Control", avg_topic_by_mouse,NA)) %>%
  mutate(control_mean = mean(ctrl_topic_values, na.rm = TRUE)) %>%
  mutate(RQ = avg_topic_by_mouse - control_mean) 

##### barplot summary - ribo removed for final figure
ggplot(subset(MITRG_plot_df, topic != "Foxn1_Ribo_topic1"), 
       aes(x = topic, y = RQ, fill = multiseq.barcodes.adjusted)) + 
  theme_bw() + geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = c("Bar1" = "cadetblue1", 
                               "Bar2" = "cadetblue3", 
                               "Bar3" = "cadetblue4", 
                               "Bar4" = "firebrick1", 
                               "Bar5" = "firebrick3", 
                               "Bar6" = "firebrick4")) + 
  ylab("Relative Foxn1 topic score") + RotatedAxis()
