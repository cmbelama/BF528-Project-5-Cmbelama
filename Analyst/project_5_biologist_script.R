
library(tidyverse)
library(dplyr)

######## 7.1 ##########

# import all fpkm tables for samples from programmer
# order by tracking id
# change fpkm column in each dataset to corresponding sample name for plotting

P0_1 <- read.table("/projectnb/bf528/users/frizzled/project_2/Programmer/528_Project_2_FINAL/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

P0_2 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/P02genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

P4_1 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/P41genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

P4_2 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/P42genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

P7_1 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/P71genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

P7_2 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/P72genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

AD_1 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/Ad1genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

AD_2 <- read.table("/projectnb/bf528/users/frizzled/project_2/Biologist/FPKM_Tables/Ad2genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)

colnames(P0_1)[10] <- "P0_1_FPKM"
colnames(P0_2)[10] <- "P0_2_FPKM"
colnames(P4_1)[10] <- "P4_1_FPKM"
colnames(P4_2)[10] <- "P4_2_FPKM"
colnames(P7_1)[10] <- "P7_1_FPKM"
colnames(P7_2)[10] <- "P7_2_FPKM"
colnames(AD_1)[10] <- "AD_1_FPKM"
colnames(AD_2)[10] <- "AD_2_FPKM"

# recreated figure 1d with fpkm tracking tables
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(gridExtra)
library(ggpubr)
library(dendextend)

# genes
sarcomere <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
mitochondria <- c("Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
cell_cycle <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27",
                "Cdc45", "Rad51", "Aurkb", "Cdc23")

# function to plot
color_pallette <- colorRampPalette(brewer.pal(8, "Dark2"))(16)
plot_fpkm <- function(genes, plot_title, ltag) {
  P0 <- (P0_1[P0_1$gene_short_name %in% genes, "P0_1_FPKM"] + P0_2[P0_2$gene_short_name %in% genes, "P0_2_FPKM"]) / 2
  P4 <- (P4_1[P4_1$gene_short_name %in% genes, "P4_1_FPKM"] + P4_2[P4_2$gene_short_name %in% genes, "P4_2_FPKM"]) / 2
  P7 <- (P7_1[P7_1$gene_short_name %in% genes, "P7_1_FPKM"] + P7_2[P7_2$gene_short_name %in% genes, "P7_2_FPKM"]) / 2
  AD <- (AD_1[AD_1$gene_short_name %in% genes, "Ad_1_FPKM"] + AD_2[AD_2$gene_short_name %in% genes, "Ad_2_FPKM"]) / 2
  all_fpkm <- cbind(P0, P4, P7, AD)  #merged across
  row.names(all_fpkm) <- genes
  all_fpkm <- data.frame(t(all_fpkm))  #transpose
  all_fpkm$timeperiod <- rownames(all_fpkm)
  all_fpkm <- melt(all_fpkm, id.vars=c("timeperiod"))
  thm = theme(plot.title = element_text(hjust=0.5))
  ggplot(all_fpkm, aes(x=factor(timeperiod, level=c("P0", "P4", "P7", "Ad")), 
                                  y=value, colour=variable, group=variable)) + geom_line() + geom_point() + labs(title = plot_title, tag=ltag, 
                                  x="Maturation Time Stages", y="FPKM", colour="Genes") + theme_bw() + thm + scale_color_manual(values=color_pallette)
}

# now plot by each function
sarcomere_plot <- plot_fpkm(sarcomere, "Sarcomere", "A")
# Mpc1 not in our FPKM tables - removed
mitochondria_plot <- plot_fpkm(mitochondria, "Mitochondria", "B")
# Bora not in our FPKM tables - removed
cell_cycle_plot <- plot_fpkm(cell_cycle, "Cell Cycle", "C")

# combine and arrange plots
combined_plot <- grid.arrange(sarcomere_plot, mitochondria_plot, cell_cycle_plot, nrow=3)
annotate_figure(combined_plot, top= "FPKM Values of Sarcomere, Mitochondira, and Cell Cycle Genes During Maturation Stages\n", fig.lab.size = 10)


######### 7.3 ############

# Merge tables
P0_1_FPKM <- P0_1[c(1, 5, 10)]
full_fpkm_matrix <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", header=TRUE, sep="\t") %>% distinct()
fpkm_combined <- merge(P0_1_FPKM, fpkm_mat, by="tracking_id")

# Identify and remove duplicated rows
duplicated <- fpkm_combined %>% group_by(tracking_id) %>% dplyr::filter(n() > 1)
`%notin%` <- Negate(`%in%`) # function to match real with duplicated
fpkm_combined = subset(fpkm_combined, fpkm_combined$tracking_id %notin% duplicated$tracking_id)
add_duplicated <- duplicated %>% group_by(tracking_id) %>% dplyr::filter(n() > 4)
duplicated <- subset(duplicated, duplicated$tracking_id %notin% add_duplicated$tracking_id)

# Extract the correct number of rows in each
duplicated <- duplicated %>% group_by(tracking_id) %>% slice(c(1,4))
add_duplicated <- add_duplicated[c(3, 5, 7, 11, 14, 20, 25), ]

# Add duplicated rows back into merged dataset
fpkm_combined <- rbind(fpkm_combined, duplicated, add_duplicated)

# load in differentially expressed genes
diff_exp <- read.table("gene_exp.diff", header = TRUE)

# find top 750 differentially expressed genes and order by q-value
diff_exp <- diff_exp %>% dplyr::filter(significant=="yes")
top_diff_exp <- head(diff_exp[order(diff_exp$q_value), ], n=750)

# Some rows have more than one gene symbol, we only want the first one
first_gene <- c(which(grepl(",", top_diff_exp$gene))) # chooses first before comma
gene_symbol <- top_diff_exp[-first_gene, c("gene")]
x <- vector()
for (i in first_gene) {
  x <- append(v, str_split(top_diff_exp[i, "gene"], ",")[[1]][1])
}

# Get the gene symbols
gene_symbol <- c(gene_symbol, x)

#Option 4 (taken)
#As the above methods returned differing numbers of Ensembl ids, I utilized the gene symbol
#Using Ensembl tracking_ids would have been an ideal method within the above Bioconductor packages, however, there were inconsistencies with the tracking_ids that I observed

# Subset this matrix by at most the top 750 genes found to be differentially expressed between P0 and Ad from 5.4
fpkm_combined_sub <- subset(fpkm_combined, fpkm_combined$gene_short_name %in% gene_symbol)

# Remove the tracking_id and gene_short_name columns and add gene names back 
final <- fpkm_combined_sub[-1:-2]
rownames(final) <- make.unique(fpkm_combined_sub[,2], sep="*") # add * when there was a duplicate found

# Some rows had 0, and needed to be fixed for clustered heatmap
final <- final[apply(final, 1, function(x) !all(x==0)), ]
final <- as.matrix(final)

# heatmap
colorss <- rev(colorRampPalette(c("purple", "black", "yellow"))(100))  #color palette
plt <- pheatmap::pheatmap(final)[[2]]
plt <- dendextend::rotate(plt, order = c("P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", 
                                         "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "Ad_1_FPKM", 
                                         "Ad_2_FPKM"))
png(filename="heatmap.png", width=950, height=950, res=200)
pheatmap::pheatmap(final, col=colorss, scale="row", cluster_rows = TRUE, 
                   clustering_distance_rows = "euclidean", cluster_cols=as.hclust(plt), 
                   show_rownames = FALSE, legend=TRUE, angle_col = 90, fontsize_col = 10)
dev.off()
