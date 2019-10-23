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
# Modify the supplied Targets.txt file by entering the file names of the .crop.txt files in column one, the group name of the corresponding sample (e.g. WT) in column 2, the time point in column 3 (if no time points were done add t0 in all rows) and the replicate number in column 4. The file must be tab separated.

sample_info <- read.delim("Targets.txt")
targets <- data.frame(files = sample_info$files, group = paste0(sample_info$group, sample_info$time),
                      description = paste0(sample_info$group, sample_info$time, "_", sample_info$replicate))
row.names(targets) <- targets$description

d <- readDGE(targets, header = FALSE) # this reads in the count data
#d$samples
#head(d$counts)
#summary(d$counts)
#dim(d)

########## Filter out non-/lowly expressed genes ########################################################
keep <- filterByExpr(d) ### from edgeR manual
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis
plotMDS(d, col=as.numeric(d$samples$group)) # plot leading fold changes - to get an idea if data have batch effects

# Uncomment the next 3 lines for exporting the MDS plot
#png(filename = "MDSplot.png", res = 300, width = 1920, height = 1440)
#plotMDS(d, col=rep(1:length(targets[,2]), each=(length(targets[,2])/length(unique(targets[,2])))))
#dev.off()

# Model fitting without batch effect correction ##########################################################################
# create design matrix
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
my.contrasts <- makeContrasts(
  t0hir1vswt = hir1t0-wtt0,
  t30hir1vswt = hir1t30-wtt30,
  wtt30vst0 = wtt30-wtt0,
  hir1t30vst0 = hir1t30-hir1t0,
  t30hir1vswtvst0 = (hir1t30-hir1t0)-(wtt30-wtt0),
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

qlf <- glmQLFTest(fit, contrast=my.contrasts[,c("t0hir1vswt")]) # differential expression analysis with the comparison (=contrast) specified in quotes (multiple comparisons can be done separated by a comma)
#qlf <- glmQLFTest(fit, contrast=c(1,-1,0,0,0,0)) ### same as above but in other format - see edgeR manual
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

# Export results
write_tsv(diff_results, "diff_results.tsv")


### C. albicans annotation data import ##########################################################
# The following section adds gene names and Entrez IDs to the result table.
library(GenomicFeatures)

gff_file <- "../required_files/C_albicans_SC5314_A22_current_features_haploid.gff"
GeneData <- rtracklayer::import.gff(gff_file)

CGDID_GeneName <- elementMetadata(GeneData)[GeneData$type == "gene", c(5,12)]
CGDID_GeneName$locus_tag <- gsub(pattern = "_", replacement = "", x = CGDID_GeneName$ID)

gff_file2 <- "../required_files/GCF_000182965.3_ASM18296v3_genomic.gff"
GeneData2 <- rtracklayer::import.gff(gff_file2)

CGDID_Entrez <- as.data.frame(elementMetadata(GeneData2[GeneData2$type == "gene", c(20,6)]))
CGDID_Entrez$locus_tag <- gsub(pattern = "CAALFM_", replacement = "", x = CGDID_Entrez$locus_tag)
CGDID_Entrez$Dbxref <- gsub(pattern = "GeneID:", replacement = "", x = CGDID_Entrez$Dbxref)

CaGeneIDs <- base::merge(CGDID_GeneName, CGDID_Entrez, by.x=3, by.y=1, all=TRUE)
colnames(CaGeneIDs) <- c("locus_tag", "CGDID", "GeneName", "EntrezID")



diff_results2 <- merge(diff_results, CaGeneIDs, by.x=0, by.y=2, all.x=TRUE)

# Export results with gene names
write_tsv(diff_results2, "diff_results_wNames.tsv")



######################## Model fitting with batch effect correction ####################################

Strain <- sample_info$group
Strain <- relevel(Strain, ref = "wt")
Time <- sample_info$time
Repl <- sample_info$replicate

design_paired <- model.matrix(~Time+Time:Repl+Time:Strain)
#design_paired <- model.matrix(~Strain+Strain:Repl+Strain:Time) # different setup - see edgeR manual
rownames(design_paired) <- colnames(d)

dp <- estimateDisp(d, design_paired, robust=TRUE) ## does both dispersion estimations in one step
plotBCV(dp) # plots the results of the dispersion estimation

# Uncomment the next 3 lines to export the dispersion plot
#pdf(file = "disp.pdf")
#plotBCV(dp)
#dev.off()

fitp <- glmQLFit(dp, design_paired, robust = TRUE) # fits a model onto the data
plotQLDisp(fitp) # plots the result of the model fitting

qlfp <- glmQLFTest(fitp, coef = c("Timet0:Strainhir1")) # performs the differential expression analysis
#qlfp <- glmQLFTest(fitp, coef = "Timet30:Strainhir1") # alternative for testing for one condition
#qlfp <- glmQLFTest(fitp, contrast = c(0,0,0,0,0,0,0,-1,1)) # to test for time point specificity
topTags(qlfp) # shows the top differentially expressed genes
summary(decideTests(qlfp)) # shows a summary of the diff. expr. analysis
plotMD(qlfp) # plots the results

# Uncomment the next 4 lines to export the results plot
#pdf("diff_expr_paired_results.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotMD(qlfp)
#abline(h=cutoff, col="blue")
#dev.off()

diff_results_paired <- as.data.frame(topTags(qlfp, n=Inf)) # saves the results in the variable diff_results_paired
diff_results_paired2 <- merge(diff_results_paired, CaGeneIDs, by.x=0, by.y=2, all.x=TRUE)

# Export results
write_tsv(diff_results_paired, "diff_results_paired.tsv")



#################### Principal component analysis ###############################################
library(ggbiplot)

cpmill <- cpm(d$counts, normalized.lib.size=TRUE)
cpmill_transp <- data.frame(t(cpmill))
colnames(cpmill_transp) <- row.names(cpmill)
cpmill_transp$group <- gsub(rownames(cpmill_transp), pattern = "_[1-3]", replacement = "")
cpmill_transp.pca <- prcomp(dplyr::select(cpmill_transp, -group), scale. = TRUE)
#PCs <- as.data.frame(cpmill_transp.pca$rotation) # exports principal components
#PCs <- PCs[order(PCs$PC1, decreasing = TRUE),] # sorts principal components

plot <- ggbiplot(cpmill_transp.pca, choices = c(1,2), var.axes = FALSE, groups = cpmill_transp$group)
plot
plot2 <- plot + cleanup + theme(legend.key = element_blank()) +
  #xlab("PC1 (XX% explained variance)") + 
  #ylab("PC2 (XX% explained variance)") + 
  scale_color_brewer(palette = "Dark2", name = expression("Strains"))
plot2  

# Uncomment these 3 lines to export plot
#pdf(file = "PCA.pdf")
#plot2
#dev.off()



#################### GO term enrichment analysis ###############################################
library(clusterProfiler)
library(AnnotationHub)
CalbicansEntrez <- AnnotationHub()[['AH67641']]


up <- diff_results2[diff_results2$FDR < 0.05 & diff_results2$logFC > 0,] # change logFC if necessary
FC <- up$logFC
names(FC) <- up$EntrezID

go <- enrichGO(up$EntrezID, OrgDb = CalbicansEntrez, readable = TRUE)

dplot <- dotplot(go, showCategory =20) + 
  scale_color_gradient(name = "Adj. p-value", low = "#cb181d", high = "#fee0d2") +
  labs(size = "No. of Genes")
dplot

cplot <- cnetplot(go, showCategory = 15,  foldChange = FC) +
  scale_color_gradient(low = "#ccece6", high = "#00441b") +
  labs(color = "Fold Change\n(log2)", size = "No. of Genes")
cplot



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
                                                  as.character(GeneName),'')),
                                 size = 3  ,hjust=-0,vjust=-0.5, colour="black")
myplot_lab

# Uncomment next 3 lines to export plot
#pdf('MyPlot.pdf')
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

vulcplot_lab <- vulcplot + geom_text(aes(label=ifelse(((logFC >= 2 & FDR < 0.05) | (logFC <= -2 & FDR < 0.05)),
                                                  as.character(GeneName),'')),
                                 size = 3  ,hjust=-0,vjust=-0.5, colour="black")
vulcplot_lab

# Uncomment next 3 lines to export plot
#pdf('MyPlot.pdf')
#vulcplot_lab
#dev.off()
