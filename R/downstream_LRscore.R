# parse input parameters from linux scripts
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Please supply the results_folder, giotto_object_path!", call.=FALSE)
} 

results_folder = args[1]
giotto_object_path = args[2]



#---------------- load packages ----------------#
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("welch-lab/cytosignal")
devtools::install_github("lhe17/nebula", force = TRUE)

library(Giotto)
library(data.table)
library(Matrix)
library(mltools)
library(cytosignal)
library(nebula)
library(ggplot2)
library(pheatmap)

runRawScore <- function(dge.path, cells.loc.path, clusters.path, prefix, res_path){
  dge = readRDS(dge.path)
  cells.loc = as.matrix(readRDS(cells.loc.path))
  colnames(cells.loc) = c("x", "y")
  clusters = readRDS(clusters.path)

  cells.loc = cells.loc[colnames(dge), ]
  clusters = clusters[colnames(dge)]
  
  # all.equal( rownames(cells.loc), colnames(dge), check.attributes=FALSE )

  cs <- createCytoSignal(raw.data = dge,
                         cells.loc = cells.loc,
                         clusters = clusters
  )
  cs <- addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)

  # remove all objects that are not in the dge
  # rm(dge, cells.loc, barcodes, genes, loc.rownames, clust.names, cell.order)

  cs <- removeLowQuality(cs, counts.thresh = 0)
  cs <- changeUniprot(cs)

  #####--------- LRscore ---------#####
  cs <- inferEpsParams(cs, scale.factor = 0.85)
  cs@parameters$r.diffuse.scale
  cs@parameters$sigma.scale

  #####--------- LRscore ---------#####
  # set the weight of index cell of diff to be 1
  cs <- findNN(cs, diff.weight = 1)
  # View(cytosignal:::findNN)
  median(table(cs@imputation[["GauEps"]]@nn.id))
  median(table(cs@imputation[["DT"]]@nn.id))

  cs <- imputeLR(cs)

  # infer LRscore without normalizing
  cs <- inferScoreLR(cs, lig.slot = "GauEps", recep.slot = "DT",
                     norm.method = "none", intr.db.name = "diff_dep")
  cs <- inferScoreLR(cs, lig.slot = "DT", recep.slot = "DT",
                     norm.method = "none", intr.db.name = "cont_dep")
  
  tmp1 = cs@lrscore[["GauEps-DT"]]@score
  print(summary(tmp1@x[tmp1@x > 1]))
  
  #####--------- rewrite spatial statistics into matrix form ---------#####
  diff.score = cs@lrscore[["GauEps-DT"]]@score
  colnames(diff.score) = unname(cytosignal:::getIntrNames(cs, colnames(diff.score)))

  cont.score = cs@lrscore[["DT-DT"]]@score
  colnames(cont.score) = unname(cytosignal:::getIntrNames(cs, colnames(cont.score)))

  # take only the parts start with "E" from rds
  # rds = substr(rds, regexpr("E", rds), nchar(rds))

  clusters = as.character(cs@clusters)
  names(clusters) = names(cs@clusters)

  saveRDS(diff.score, paste0(res_path, "/DT-Diff_LRscore-", prefix, ".RDS"))
  saveRDS(cont.score, paste0(res_path, "/DT-Cont_LRscore-", prefix, ".RDS"))
  saveRDS(clusters, paste0(res_path, "/Clusters-", prefix, ".RDS"))

  gc()
}

concatMetaData <- function(data.df, all.cluster, rds_list){

  res.list = lapply(seq_along(data.df$sid), function(x){
    sid = data.df$sid[x]
    # rds_path = rds_list[grepl(paste0(fov, "-RAW-", cat), rds_list)]
    scores = readRDS(rds_list[x])
    # for count that larger than 0 but less than 1, set it to 0
    # for counts that larger than 1, round them to integer
    scores@x[scores@x < 1] <- 0
    scores@x <- round(scores@x)
    # re-format scores as dgCMatrix
    scores = as(scores, "matrix")
    scores = as(scores, "dgCMatrix")

    rownames(scores) = paste0("d_", x, "-", rownames(scores))

    # get the cell type
    # cell_type = all.cluster[grepl(fov, names(all.cluster)), drop = T]
    # cell_type = cell_type[rownames(scores)]
    # cell_type = as.character(cell_type)
    cell_type = all.cluster[rownames(scores)]

    # get the metadata
    class_vec = rep(data.df$class[data.df$sid == sid], nrow(scores))

    # get the rep
    sid_vec = rep(data.df$sid[data.df$sid == sid], nrow(scores))

    return.list = list(
      scores = scores,
      cell_type = cell_type,
      class = class_vec,
      sid = sid_vec)
    names(return.list) = c("scores", "cell_type", "class", "sid")

    return(return.list)
  })

  # iterate over the scores and intersect colnames
  intr.list = Reduce(intersect, lapply(res.list, function(x) colnames(x$scores)))
  all.scores = lapply(res.list, function(x){
    x$scores = x$scores[,intr.list]
    return(x$scores)
  })
  all.scores = Reduce(rbind, all.scores)
  all.cell_type = unlist(lapply(res.list, function(x) x$cell_type))
  all.class = unlist(lapply(res.list, function(x) x$class))
  all.sid = unlist(lapply(res.list, function(x) x$sid))

  pred = data.frame(
    cell_type = all.cell_type,
    class = all.class
  )

  return(list(
    count = t(as.matrix(all.scores)),
    pred = pred,
    sid = all.sid
  ))
}

runNEBULA <- function(
    diff.data,
    cont.data,
    diff.model,
    cont.model,
    save.dir
){

  #---------------------------------------- run diff ----------------------------------------#
  cat("Running Diffusion-dependent interactions...\n")
  # running nebula
  diff.re <- nebula(diff.data$count, diff.data$sid, pred=diff.model, ncore=8)
  print(paste0("Non-Complete cases:", nrow(diff.re$summary) - sum(complete.cases(diff.re$summary))))

  # remove rows with NA
  diff.summary = diff.re$summary[complete.cases(diff.re$summary),]

  # For each predictor, output the top 5, middle 2, and bottom 5 gene, as a dataframe
  intrs.use = colnames(diff.summary)[grepl("p_", colnames(diff.summary))][-1]
  diff.res.rank = sapply(intrs.use, function(x){
    tmp = diff.summary[order(diff.summary[,x]), "gene"]
    top5 = tmp[1:5]
    mid2 = tmp[(length(tmp)/2-1):(length(tmp)/2)]
    bot5 = tmp[(length(tmp)-4):length(tmp)]
    return(data.frame(c(top5, mid2, bot5)))
  }, USE.NAMES=T, simplify = F)

  diff.res.rank = Reduce(cbind, diff.res.rank)
  colnames(diff.res.rank) = intrs.use

  saveRDS(diff.summary, file.path(save.dir, "NEBULA_diff_summary.rds"))
  # saveRDS(diff.res.rank, file.path(save.dir, "NEBULA_diff_rank.rds"))
  # saveRDS(diff.re, file.path(save.dir, "NEBULA_diff_all_results.rds"))

  # how many genes are predicted sigifinciant for each predictor
  diff.pval = diff.summary[, c("gene", grep("p_", colnames(diff.summary), value=T) )]
  rownames(diff.pval) = diff.pval$gene

  if(ncol(diff.pval) > 3){
    diff.pval = diff.pval[, c(-1, -2)]
  } else{
    g.names = diff.pval$gene
    p.name = colnames(diff.pval)[length(colnames(diff.pval))]
    diff.pval = diff.pval[, c(-1, -2)]
    diff.pval = data.frame(
      diff.pval
    )
    colnames(diff.pval) = p.name
    rownames(diff.pval) = g.names
  }

  # apply FDR correction for each column
  diff.pval = apply(diff.pval, 2, function(x){
    p.adjust(x, method = "fdr")
  })

  diff.prop = apply(diff.pval < 0.05, 2, sum) / nrow(diff.pval)
  diff.prop = as.data.frame(diff.prop)
  colnames(diff.prop) = "significant_proportion"
  # remove "p_celltype" prefix
  rownames(diff.prop) = gsub("p_cell_type", "", rownames(diff.prop))

  p1 = ggplot(data.frame(diff.prop),
              aes(x=reorder(rownames(diff.prop),-significant_proportion),
                  y=significant_proportion)) +
    # sort the barplot by the proportion of significant genes
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title=paste0(
      "Proportion of significant LR for each predictor (perturb% = ",
      round(diff.prop["p_classPerturb", ], 3), ")"),
      x="Predictor", y="Proportion of significant genes")
  # p1

  #---------------------------------------- run cont ----------------------------------------#
  cat("Running Contact-dependent interactions...\n")
  cont.re <- nebula(cont.data$count, cont.data$sid, pred=cont.model, ncore=5)
  print(paste0("Complete cases:", nrow(cont.re$summary) - sum(complete.cases(cont.re$summary))))

  cont.summary = cont.re$summary[complete.cases(cont.re$summary),]

  intrs.use = colnames(cont.summary)[grepl("p_", colnames(cont.summary))][-1]
  cont.res.rank = sapply(intrs.use, function(x){
    tmp = cont.summary[order(cont.summary[,x]), "gene"]
    top5 = tmp[1:5]
    mid2 = tmp[(length(tmp)/2-1):(length(tmp)/2)]
    bot5 = tmp[(length(tmp)-4):length(tmp)]
    return(data.frame(c(top5, mid2, bot5)))
  }, USE.NAMES=T, simplify = F)

  cont.res.rank = Reduce(cbind, cont.res.rank)
  colnames(cont.res.rank) = intrs.use

  saveRDS(cont.summary, file.path(save.dir, "NEBULA_cont_summary.rds"))
  # saveRDS(cont.res.rank, file.path(save.dir, "NEBULA_cont_rank.rds"))
  # saveRDS(cont.re, file.path(save.dir, "NEBULA_cont_all_results.rds"))

  cont.pval = cont.summary[, c("gene", grep("p_", colnames(cont.summary), value=T) )]
  rownames(cont.pval) = cont.pval$gene

  if(ncol(cont.pval) > 3){
    cont.pval = cont.pval[, c(-1, -2)]
  } else{
    g.names = cont.pval$gene
    p.name = colnames(cont.pval)[length(colnames(cont.pval))]
    cont.pval = cont.pval[, c(-1, -2)]
    cont.pval = data.frame(
      cont.pval
    )
    colnames(cont.pval) = p.name
    rownames(cont.pval) = g.names
  }

  cont.pval = apply(cont.pval, 2, function(x){
    p.adjust(x, method = "fdr")
  })

  cont.prop = apply(cont.pval < 0.05, 2, sum) / nrow(cont.pval)
  cont.prop = as.data.frame(cont.prop)
  colnames(cont.prop) = "significant_proportion"
  rownames(cont.prop) = gsub("p_cell_type", "", rownames(cont.prop))

  p2 = ggplot(data.frame(cont.prop),
              aes(x=reorder(rownames(cont.prop), -significant_proportion),
                  y=significant_proportion)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(
      "Proportion of significant LR for each predictor (perturb% = ",
      round(cont.prop["p_classPerturb", ], 3), ")"),
      x="Predictor", y="Proportion of significant genes")

  png(file.path(save.dir, "ALL_prop.png"), width=8, height = 8, units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, nrow=2))
  dev.off()

  cat("Done!\n")

}



#---------------- Prepare Data for CytoSignal ----------------#
# set working directory
setwd(results_folder)

# psoriasis
xenium_gobj = loadGiotto(giotto_object_path)
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)

xenium_gobj = replaceGiottoInstructions(xenium_gobj, instructions = instrs)

cell_info = pDataDT(xenium_gobj)[, c('cell_ID', 'list_ID', 'leiden_clus', 'cell_type_isML_subcelltyping')]

# create dge file as tutorial
raw_mat = as.matrix(xenium_gobj@expression$cell$rna$raw@exprMat)
raw_mat = t(raw_mat)
x <- data.table(raw_mat)
sparseM <- sparsify(x)
rownames(sparseM) <- rownames(raw_mat) 
head(sparseM)
print(summary(colSums(as.matrix(t(raw_mat)))))
##saveRDS(sparseM,'raw_count_mat.rds')

# create spatial file as tutorial
pos_file = xenium_gobj@spatial_locs$cell$raw@coordinates
pos_file = as.data.frame(pos_file)
rownames(pos_file) = pos_file$cell_ID
pos_file$cell_ID=NULL
pos_file=as.matrix(pos_file)
head(pos_file)
str(pos_file)
saveRDS(pos_file,'count_spatial.rds')

# create cluster file as tutorial
celltype = factor(cell_info[(cell_info$list_ID=="obj1")|(cell_info$list_ID=="obj2")|(cell_info$list_ID=="obj3")|(cell_info$list_ID=="obj4"),]$cell_type_isML_subcelltyping)
names(celltype) = rownames(raw_mat)
str(celltype)
saveRDS(celltype,'clusters_nebula.rds')

# select control and trt group
head(a1_cell_ID,n=20)
a1_cell_ID <- t(sparseM[grepl("^obj1", rownames(sparseM)), ])
a2_cell_ID <- t(sparseM[grepl("^obj2", rownames(sparseM)), ])
b_cell_ID <- t(sparseM[grepl("^obj3", rownames(sparseM)), ])
c1_cell_ID <- t(sparseM[grepl("^obj4", rownames(sparseM)), ])
saveRDS(a1_cell_ID,'a1_count_mat.rds')
saveRDS(a2_cell_ID,'a2_count_mat.rds')
saveRDS(b_cell_ID,'b_count_mat.rds')
saveRDS(c1_cell_ID,'c1_count_mat.rds')



#---------------- Getting raw LRscore without normalization ----------------#
# run this step for all datasets you have
pos_files <- c("a1", "a2", "b", "c1")

# Create output path
res_path = file.path(results_folder, "multi-sample")
if (!dir.exists(res_path)) {
  dir.create(res_path)
}

# Loop through pos_files
for (pos in pos_files) {
  print(paste("Processing:", pos, "data file..."))
  dge.path = file.path(results_folder, paste0(pos, "_count_mat.rds")) # gene X cell matrix path
  
  # Run CytoSignal, testing removeLowQuality and inferEpsParams settings
  runRawScore(
    dge.path,
    file.path(results_folder, "count_spatial.rds"),  # Cell location path
    file.path(results_folder, "clusters_nebula.rds"),  # Cluster assignment path
    pos,  # Current pos ("a1", "a2", "b", "c1")
    res_path  # Output path
  )
  
  # Print message after each file is processed
  print(paste(pos, "data file processing complete!"))
}


#---------------- Run NEBULA ----------------#
# Need to change the settings based on how many datasets you have for each
# condition, and how many replicates you have for each dataset

# List of prefixes for different datasets
prefixes_control = c("a1", "a2")
prefixes_treatment = c("b", "c1")
# Define control paths dynamically based on prefix
control.diff.path1 = file.path(res_path, paste0("DT-Diff_LRscore-", prefixes_control[1], ".RDS"))
control.cont.path1 = file.path(res_path, paste0("DT-Cont_LRscore-", prefixes_control[1], ".RDS"))
control.cluster.path1 = file.path(res_path, paste0("Clusters-", prefixes_control[1], ".RDS"))

control.diff.path2 = file.path(res_path, paste0("DT-Diff_LRscore-", prefixes_control[2], ".RDS"))
control.cont.path2 = file.path(res_path, paste0("DT-Cont_LRscore-", prefixes_control[2], ".RDS"))
control.cluster.path2 = file.path(res_path, paste0("Clusters-", prefixes_control[2], ".RDS"))

# Define treatment paths dynamically based on prefix
treatment.diff.path1 = file.path(res_path, paste0("DT-Diff_LRscore-", prefixes_treatment[1], ".RDS"))
treatment.cont.path1 = file.path(res_path, paste0("DT-Cont_LRscore-", prefixes_treatment[1], ".RDS"))
treatment.cluster.path1 = file.path(res_path, paste0("Clusters-", prefixes_treatment[1], ".RDS"))

treatment.diff.path2 = file.path(res_path, paste0("DT-Diff_LRscore-", prefixes_treatment[2], ".RDS"))
treatment.cont.path2 = file.path(res_path, paste0("DT-Cont_LRscore-", prefixes_treatment[2], ".RDS"))
treatment.cluster.path2 = file.path(res_path, paste0("Clusters-", prefixes_treatment[2], ".RDS"))

# Ensure the output directory exists
save.dir = file.path(res_path, "nebula_results")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}


diff_rds_list = c(
  control.diff.path1,control.diff.path2,
  treatment.diff.path1,treatment.diff.path2
)

cont_rds_list = c(
  control.cont.path1,control.cont.path2,
  treatment.cont.path1,treatment.cont.path2
)

clusters.path = c(
  control.cluster.path1,control.cluster.path2,
  treatment.cluster.path1,treatment.cluster.path2
)

clusters.list = lapply(seq_along(clusters.path), function(x){
  tmp1 = readRDS(clusters.path[x])
  names(tmp1) = paste0("d_", x, "-", names(tmp1))
  return(tmp1)
})
all.cluster = Reduce(c, clusters.list)

meta.df = data.frame(
  sid = c("1", "2","3","4"),
  class = c("Ctrl", "Ctrl","Treatment","Treatment")
)

diff.data = concatMetaData(meta.df, all.cluster, diff_rds_list)
cont.data = concatMetaData(meta.df, all.cluster, cont_rds_list)

# predictor: cell_type + class
diff.model = model.matrix(~cell_type + class, data=diff.data$pred)
cont.model = model.matrix(~cell_type + class, data=cont.data$pred)

# predictor: class
diff.model = model.matrix(~class, data=diff.data$pred)
cont.model = model.matrix(~class, data=cont.data$pred)

runNEBULA(
  diff.data,
  cont.data,
  diff.model,
  cont.model,
  save.dir
)


#---------------- Distinguish if the significant LRscore are higher or lower in cases/controls----------------#
# Load the summary data using the save.dir in NEBULA
summaries = list(
  "diff" = file.path(save.dir, "NEBULA_diff_summary.rds"),
  "cont" = file.path(save.dir, "NEBULA_cont_summary.rds")
)

# Define plot titles and file names for each summary type
plot_titles = list(
  "diff" = "Volcano Plot for diffuse interaction",
  "cont" = "Volcano Plot for contact interaction"
)

# Start the loop to process each summary type
for (type in names(summaries)) {
  # Load the summary data
  summary_data=  readRDS(summaries[[type]])

  # Define the output path for the PNG file
  output_path = file.path(save.dir, paste0("volcano_plot_", type, ".png"))
  
  # Create the volcano plot
  png(output_path, width = 800, height = 600)
  plot(summary_data$logFC_classTreatment, -log10(summary_data$p_classTreatment),
       xlab = "logFC", ylab = "-log10(p-value)",
       main = plot_titles[[type]], pch = 19, col = "blue")
  
  # Add significance threshold line
  abline(h = -log10(0.05), col = "red", lty = 2)
  
  # Identify significant points based on p-value and logFC
  significant <- summary_data$p_classTreatment < 0.05 & abs(summary_data$logFC_classTreatment) > 1
  
  # Label significant points with gene names (or other appropriate labels)
  text(summary_data$logFC_classTreatment[significant], 
       -log10(summary_data$p_classTreatment[significant]), 
       labels = summary_data$gene[significant],  # Change "gene" to actual column name if needed
       pos = 4, cex = 0.7, col = "black")
  
  # Close the PNG device
  dev.off()
}


#---------------- cell types contribution ----------------#
## ALL_prop plot shows the proportion of LR with p value less than 0.05 for each cell types

