# Project 5 - Analyst - Camilla Belamarich

# Load in Differential Expression sample data - P4 vs. P7
diff_exp <- read.table("gene_exp.diff", header = TRUE)

# Sort the data frame so that the smallest q_values are at the top
sorted_diff_exp <- diff_exp[order(diff_exp[,13]),]

# Produce a table of the top ten differentially expressed genes, with their names, 
# FPKM values, log fold change, p-value, and q-value
top_ten_genes <- sorted_diff_exp[1: 10,]
columns_interested <- top_ten_genes[,c(3,8,9,10,12,13)] # only the columns we are interested in
# rename columns
names(columns_interested) <- c('Gene Name','FPKM Value_1','FPKM Value_2','log_fold change','P-Value','Q-Value')

# write table of top 10 differentiall expressed genes
write.csv(columns_interested,"top_10_de_genes.csv")

# Produce a histogram of the log2.foldchange column for all genes
hist(diff_exp$log2.fold_change., breaks = 40, xlab = "Log2 Fold Change",
     main = "Histogram of Log2 Fold Change for All Genes")
nrow(diff_exp) # 36329

# Create a new data frame that contains only the significant genes
sig_genes <- subset(diff_exp, diff_exp$significant == "yes")
nrow(sig_genes) # 2139 total significant genes
sig_up2 <- subset(sig_genes,log2.fold_change.>0)
sig_down2 <- subset(sig_genes,log2.fold_change.<0)
nrow(sig_up2) # 1084
nrow(sig_down2) # 1055

# Create a second histogram of the log2 fold change values only for significant genes
hist(sig_genes$log2.fold_change., breaks = 30, xlab = "Log2 Fold Change",
     main = "Histogram of Log2 Fold Change for All Significant Genes")
nrow(diff_exp) # 36329
# There is no data at value 0, this shows the up- and down-regulated genes

# Further subset the significant gene data frame you just created into two
# separate data frames with only the up- and down-regulated genes using the log2.foldchange column
# p-val < 0.01
new_sig <- subset(diff_exp, diff_exp$p_value < 0.01)
nrow(new_sig) # 2376
sig_up <- subset(new_sig,log2.fold_change.>0)
sig_down <- subset(new_sig,log2.fold_change.<0)
nrow(sig_up) # 1187 sig up-regulated genes
nrow(sig_down) # 1189 down-regulated genes

# write table for up- and down- regulated genes
write.csv(sig_up,"sig_up.csv")
write.csv(sig_down,"sig_down.csv")
