######## Import libraries and plot cleanup script ##############################################
library(edgeR)
library(tidyverse)
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/plot_cleanup.R", ssl.verifypeer = FALSE)))

######### Set parameter for analysis #################################################################
pval <- 0.05 # set pvalue threshold (only relevant for plotting)
pval_adjust <- "BH"  # set method for p-value adjustment for multiple testing (= FDR calculation)
cutoff <- c(-1,1)  # set log2 fold change line in plots (e.g. "c(-1,1)" means lines will be drawn at 2-fold up and down)

########## Import count data #######################################################################
sample_info <- read.delim("Targets.txt")
targets <- data.frame(files = sample_info$files, group = paste0(sample_info$group, sample_info$time),
                      description = paste0(sample_info$group, sample_info$time, "_", sample_info$replicate))
row.names(targets) <- targets$description

d <- readDGE(targets, header = FALSE) # this reads in the count data
#d$samples
head(d$counts)
#summary(d$counts)
#dim(d)

########## Filter out non-/lowly expressed genes ########################################################
keep <- filterByExpr(d) ### from edgeR manual
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis
plotMDS(d, col=sample_info$replicate) # plot leading fold changes - to get an idea if data have batch effects

# Uncomment the next 3 lines for exporting the MDS plot
#png(filename = "MDSplot.png", res = 300, width = 1920, height = 1440)
#plotMDS(d, col=rep(1:length(targets[,2]), each=(length(targets[,2])/length(unique(targets[,2])))))
#dev.off()

# Model fitting without batch effect correction ##########################################################################
# create design matrix
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
my.contrasts <- makeContrasts(
  t0vst2 = DC_Cgt2-DCt0,
  t0vst6 = DC_Cgt6-DCt0,
  # include more comparisons here if necessary
  levels=design)

# Estimating the dispersions and plot them
d <- estimateDisp(d, design, robust=TRUE) ## does both dispersion estimations in one step - suggested in edgeR manual;

plotBCV(d) # plots the results of the dispersion estimation

# Uncomment the next 3 lines to export the dispersion plot
#pdf("disp.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotBCV(d)
#dev.off()# closes and saves the pdf file with the plot

fit <- glmQLFit(d, design, robust = TRUE) # fits a model onto the data
plotQLDisp(fit) # plots the result of the model fitting

qlf <- glmQLFTest(fit, contrast=my.contrasts[,c("t0vst6")]) # differential expression analysis with the comparison (=contrast) specified in quotes (multiple comparisons can be done separated by a comma)
#qlf <- glmQLFTest(fit, contrast=my.contrasts[,c("t0vst2","t0vst6")]) # diff expr. genes at either t2 or t6

topTags(qlf) # shows the top differentially expressed genes
summary(decideTests(qlf)) # gives a summary of the differential expression analysis

plotMD(qlf) # plots the results
abline(h=c(-1, 1), col="blue")

# Uncomment the next 4 lines to export the results plot
#pdf("diff_expr_results.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotMD(qlf)
#abline(h=cutoff, col="blue")
#dev.off()

diff_results <- as.data.frame(topTags(qlf, n=Inf)) # saves the results in the variable diff_results

library(clusterProfiler)
library(org.Mm.eg.db)

geneIDs <- bitr(rownames(diff_results), fromType = "ENSEMBL", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = org.Mm.eg.db, drop = FALSE)

diff_results2 <- merge(diff_results, geneIDs, by.x=0, by.y=1)
names(diff_results2)[1] <- "EnsemblID"

# Export results
write_tsv(diff_results2, "diff_results.tsv")



######################## Model fitting with batch effect correction ####################################

Group <- d$samples$group
Group <- relevel(Group, ref = "DCt0")
#Time <- sample_info$time
Repl <- as.factor(sample_info$replicate)

design_paired <- model.matrix(~Repl+Group) # see edgeR manual for design setup
rownames(design_paired) <- colnames(d)
colnames(design_paired) <-gsub(pattern = "Group", replacement = "", x = colnames(design_paired))

dp <- estimateDisp(d, design = design_paired, robust=TRUE) ## does both dispersion estimations in one step
plotBCV(dp) # plots the results of the dispersion estimation

# Uncomment the next 3 lines to export the dispersion plot
#pdf(file = "disp.pdf")
#plotBCV(dp)
#dev.off()

fitp <- glmQLFit(dp, design_paired, robust = TRUE) # fits a model onto the data
plotQLDisp(fitp) # plots the result of the model fitting

qlfp <- glmQLFTest(fitp, coef = "DC_Cgt6") # performs the differential expression analysis
#qlfp <- glmQLFTest(fitp, coef = "Timet30:Strainhir1") # alternative for testing for one condition
topTags(qlfp) # shows the top differentially expressed genes
summary(decideTests(qlfp)) # shows a summary of the diff. expr. analysis
plotMD(qlfp) # plots the results

# Uncomment the next 4 lines to export the results plot
#pdf("diff_expr_paired_results.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotMD(qlfp)
#abline(h=cutoff, col="blue")
#dev.off()

diff_results_paired <- as.data.frame(topTags(qlfp, n=Inf)) # saves the results in the variable diff_results_paired

diff_results_paired2 <- merge(diff_results_paired, geneIDs, by.x=0, by.y=1)
names(diff_results_paired2)[1] <- "EnsemblID"

# Export results
write_tsv(diff_results_paired2, "diff_results_paired.tsv")



#################### Principal component analysis ###############################################

library(ggbiplot)

cpmill <- cpm(d$counts, normalized.lib.size=TRUE)
cpmill_transp <- data.frame(t(cpmill))
colnames(cpmill_transp) <- row.names(cpmill)
cpmill_transp$group <- gsub(rownames(cpmill_transp), pattern = "_[1-3]", replacement = "")
cpmill_transp.pca <- prcomp(dplyr::select(cpmill_transp, -group), scale. = TRUE)
PCs <- as.data.frame(cpmill_transp.pca$rotation)
PCs <- PCs[order(PCs$PC1, decreasing = TRUE),]

plot <- ggbiplot(cpmill_transp.pca, choices = c(1,2), var.axes = FALSE, groups = cpmill_transp$group)
plot
plot2 <- plot + cleanup + theme(legend.key = element_blank()) +
  #xlab("PC1 (XX% explained variance)") + 
  #ylab("PC2 (XX% explained variance)") + 
  scale_color_brewer(palette = "Dark2", name = expression("Groups"))
plot2  

# Uncomment these 3 lines to export plot
#pdf(file = "PCA.pdf")
#plot2
#dev.off()



#################### GO term enrichment analysis ###############################################
library(clusterProfiler)
library(org.Mm.eg.db)


up <- diff_results[diff_results$logFC > 1 & diff_results$FDR < 0.05,]

go <- enrichGO(rownames(up), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", readable = TRUE)

dplot <- dotplot(go, showCategory =20) + 
  scale_color_gradient(name = "Adj. p-value", low = "#cb181d", high = "#fee0d2") +
  labs(size = "No. of Genes")
dplot

FC <- up$logFC
names(FC) <- rownames(up)
cplot <- cnetplot(go, showCategory = 15,  foldChange = FC) +
  scale_color_gradient(low = "#ccece6", high = "#00441b") +
  labs(color = "Fold Change\n(log2)", size = "No. of Genes")
cplot

# Uncomment these 3 lines to export plot
#pdf(file = "GO_plots.pdf")
#dplot
#cplot
#dev.off()


############## Scatterplot ###################################################################

### Enter axis labels ###
xlb = "Average counts per million reads (log2)"
ylb = "Fold change (log2)"

### Plotting
plot_data <- diff_results2
plot_data$Group <- ifelse((plot_data$FDR < 0.05), 'diff_reg', 'not_reg')

myplot <- ggplot(data=plot_data, aes(x=logCPM, y=logFC)) +
  geom_point() +
  geom_point(data = plot_data[plot_data$Group == 'diff_reg',], color = 'red') +
  geom_hline(yintercept = c(-1,1), color='blue', linetype="dashed") +
  #xlim(c(-5,15)) +
  #ylim(c(-10,15)) +
  xlab(xlb) +
  ylab(ylb) +
  cleanup
myplot

myplot_lab <- myplot + geom_text(aes(label=ifelse(((logFC >= 5 & FDR < 0.05) | 
                                                     (logFC <= -2.5 & FDR < 0.05)),
                                                  as.character(SYMBOL),'')),
                                 size = 3  ,hjust=-0,vjust=-0.5, colour="black")
myplot_lab

# Uncomment next 3 lines to export plot
#pdf('ScatterPlot.pdf')
#myplot_lab
#dev.off()



############## Vulcano plot #################################################################

xlb_vulc = "log2 fold change hir1 vs wt"
ylb_vulc = "-log10(adjusted p-value)"

vulcplot <- ggplot(data=plot_data, aes(x=logFC, y=-log10(FDR))) +
  geom_point() +
  geom_point(data = plot_data[plot_data$Group == 'diff_reg',], color = 'red') +
  geom_vline(xintercept = c(-1,1), color='blue', linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), color='blue', linetype="dashed") +
  #xlim(c(-5,15)) +
  #ylim(c(0,130)) +
  xlab(xlb_vulc) +
  ylab(ylb_vulc) +
  cleanup
vulcplot

vulcplot_lab <- vulcplot + geom_text(aes(label=ifelse(((logFC >= 4 & FDR < 0.05) | (logFC <= -4 & FDR < 0.05)),
                                                      as.character(SYMBOL),'')),
                                     size = 3  ,hjust=-0,vjust=-0.5, colour="black")
vulcplot_lab

# Uncomment next 3 lines to export plot
#pdf('VulcPlot.pdf')
#vulcplot_lab
#dev.off()
