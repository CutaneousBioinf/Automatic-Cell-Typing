### symphony for annotating xenium data given scRNA reference's main/sub-types
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: config")
  quit(status = 1)
}

source(args[1])

###################################################
library(Matrix, lib.loc="/home/alextsoi/R/R-4.4/lib/")
library("Seurat", lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('symphony', lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('tibble')
library('dplyr')
library("irlba")


ref_metadata <- ref_metadata %>% 
                filter(!(get(maintype_col_name) %in% remove_maintype)) %>% 
                filter(!(get(subtype_col_name) %in% remove_subtype))
ref_exp <- ref_exp[,rownames(ref_metadata)]

ref_metadata_for_sub <- ref_metadata
ref_exp_for_sub <- ref_exp
########## Build the Reference for Main Celltypes ##########
if (skip_build_ref_main) {print('Skip building reference for main celltypes.')} else {
print('Start building reference for main celltypes...')
for (p in c(save_main_ref_dir, save_main_uwot_dir, save_sub_ref_dir, save_sub_uwot_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}
# Preprocess
print('Preprocessing...')
# remove missing celltypes and downsample
if (downsample) {
    ref_metadata <- ref_metadata %>% 
                    #filter(get(subtype_col_name) != '') %>% 
                    rownames_to_column() %>% 
                    group_by(get(subtype_col_name)) %>% 
                    sample_frac(downsample_to/dim(ref_metadata)[1]) %>% 
                    column_to_rownames()
} else {
    ref_metadata <- ref_metadata %>% 
                    #filter(get(subtype_col_name) != '') %>% 
                    rownames_to_column() %>% 
                    column_to_rownames()
}
ref_exp <- ref_exp[,rownames(ref_metadata)]
# remove unexpressed genes, otherwise they may cause errors
ref_exp <- ref_exp[rowSums(ref_exp)>0,]

# Calculate and save the mean and standard deviations for each gene
print('Calculating mean and sd...')
genes_means_sds = tibble(symbol = rownames(ref_exp), mean = Matrix::rowMeans(ref_exp))
genes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp, genes_means_sds$mean)

# Scale data using calculated gene means and standard deviations
print('Scaling data...')
ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, genes_means_sds$mean, genes_means_sds$stddev, 1)

# Run SVD, save gene loadings (s$u)
print('Running SVD...')
s = irlba(ref_exp_scaled, nv = num_PC, fastpath=FALSE)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

# Run Harmony integration
ref_harmObj = harmony::HarmonyMatrix(
        data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
        meta_data = ref_metadata, ## dataframe with cell labels
        theta = theta,             ## cluster diversity enforcement
        vars_use = vars_use,    ## variable to integrate out
        nclust = 100,             ## number of clusters in Harmony model
        max.iter.harmony = 20,
        return_object = TRUE,     ## return the full Harmony model object
        do_pca = FALSE,            ## don't recompute PCs
        npcs = num_PC
)



# Compress a Harmony object into a Symphony reference
if (substring(save_main_uwot_dir,1,1) != '/') {save_main_uwot_dir = paste(getwd(),'/',save_main_uwot_dir, sep='')}
reference = symphony::buildReferenceFromHarmonyObj(
                           ref_harmObj,            # output object from HarmonyMatrix()
                           ref_metadata,           # reference cell metadata
                           genes_means_sds,     # gene names, means, and std devs for scaling
                           loadings,               # genes x PCs matrix
                           verbose = TRUE,         # verbose output
                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
                           save_uwot_path = paste(save_main_uwot_dir,'/uwot_main', sep=''))

# specify which normalization method was used to build the reference 
reference$normalization_method = 'log(CP10k+1)'

# Save Symphony reference
saveRDS(reference, paste(save_main_ref_dir, '/ref_main.rds', sep=''))
}


########## Build the Reference for Subtypes ##########
if (skip_build_ref_sub) {print('Skip building reference for sub celltypes.')} else {
print('Start building reference for sub celltypes...')
print('Counts of each main cell type:')
print(table(ref_metadata[,maintype_col_name]))

for (main_type in names(table(ref_metadata[,maintype_col_name]))){
	print(paste('Start building reference for', main_type, '...'))
    ref_metadata_sub <- ref_metadata_for_sub %>%    # notice: here is not the processed data in last step
                        #filter(get(subtype_col_name) != '') %>% 
                        filter(get(maintype_col_name) == main_type)
	if (length(table(ref_metadata_sub[,subtype_col_name])) <= 1){
		print('No subtype found. Skip.')
		next
	}
    if (downsample) {
    if (dim(ref_metadata_sub)[1]>downsample_to){
        ref_metadata_sub <- ref_metadata_sub %>% 
                        rownames_to_column() %>% 
                        group_by(get(subtype_col_name)) %>% 
                        sample_frac(downsample_to/dim(ref_metadata_sub)[1]) %>% 
                        column_to_rownames()
    }}
	ref_exp_sub <- ref_exp_for_sub[,rownames(ref_metadata_sub)]
	ref_exp_sub <- ref_exp_sub[rowSums(ref_exp_sub)>0,]

	# Calculate and save the mean and standard deviations for each gene
	print('Calculating mean and sd...')
	genes_means_sds = tibble(symbol = rownames(ref_exp_sub), mean = Matrix::rowMeans(ref_exp_sub))
	genes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp_sub, genes_means_sds$mean)

	# Scale data using calculated gene means and standard deviations
	print('Scaling data...')
	ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp_sub, genes_means_sds$mean, genes_means_sds$stddev, 1)
	print(str(ref_exp_scaled))

	# Run SVD, save gene loadings (s$u)
	print('Running SVD...')
	s = irlba(ref_exp_scaled, nv = num_PC, fastpath=FALSE)
	Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
	loadings = s$u

	# Run Harmony integration
	ref_harmObj = harmony::HarmonyMatrix(
	        data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
	        meta_data = ref_metadata_sub, ## dataframe with cell labels
	        theta = theta,             ## cluster diversity enforcement
	        vars_use = vars_use,    ## variable to integrate out
	        nclust = 100,             ## number of clusters in Harmony model
	        max.iter.harmony = 20,
	        return_object = TRUE,     ## return the full Harmony model object
	        do_pca = FALSE,            ## don't recompute PCs
            npcs = num_PC
	)

	# Compress a Harmony object into a Symphony reference
    if (substring(save_sub_uwot_dir,1,1) != '/') {save_sub_uwot_dir = paste(getwd(),'/',save_sub_uwot_dir, sep='')}
	reference = symphony::buildReferenceFromHarmonyObj(
	                           ref_harmObj,            # output object from HarmonyMatrix()
	                           ref_metadata_sub,           # reference cell metadata
	                           genes_means_sds,     # gene names, means, and std devs for scaling
	                           loadings,               # genes x PCs matrix
	                           verbose = TRUE,         # verbose output
	                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
	                           save_uwot_path = paste(save_sub_uwot_dir,'/uwot_sub_',gsub(' ','-',main_type),sep='')
	                           )

	# specify which normalization method
	reference$normalization_method = 'log(CP10k+1)'

	# Save Symphony reference
	saveRDS(reference, paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds',sep=''))
}
}


########## Set up Functions for plot ##########
library(ggplot2, lib.loc="/home/alextsoi/R/R-4.4/lib/")
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(Polychrome)
library(cowplot, lib.loc="/home/alextsoi/R/R-4.4/lib/")

plotBasic <- function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                      title = 'Query',         # Plot title
                      color.by = 'cell_type',  # metadata column name for coloring
                      facet.by = NULL,         # (optional) metadata column name for faceting
                      color.mapping = NULL,    # custom color mapping
                      legend.position = 'right', # Show cell type legend
					  save_path = './plot.png',
                      size = c(8,10)) {  
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)

    if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    
    ggsave(save_path, plot=p, width=size[2], height=size[1])
}

# to make colors more distinct
# for main celltypes
main.names <- names(table(ref_metadata[maintype_col_name]))
color.mapping.main <- createPalette(length(main.names),  c("#ff0000", "#00ff00", "#0000ff"))
names(color.mapping.main) <- main.names
# for sub celltypes
sub.names <- names(table(ref_metadata[subtype_col_name]))
color.mapping.sub <- createPalette(length(sub.names),  c("#ff0000", "#00ff00", "#0000ff"))
names(color.mapping.sub) <- sub.names


#####################################################################
###     START MAPPING TARGET DATA
print("start mapping query data...")
########## Predict Main Celltype ##########
for (p in c(output_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}
print('Start working on main celltypes.')
## load reference for main celltype
reference_main <- readRDS(paste(save_main_ref_dir, '/ref_main.rds', sep=''))
# If it cannot automatically find uwot model, manually set the path.
#reference_main$save_uwot_path <- ''

# plot distribution of celltypes in reference space
plotBasic(cbind(reference_main$meta_data, reference_main$umap$embedding), 
          title = 'Reference Cells in Reference UMAP Space of Main celltypes', 
          color.by = maintype_col_name,
          color.mapping = color.mapping.main,
          save_path = paste(output_dir, '/1-main_celltype_refer.png', sep=''),
          size = c(8,10))
#plotBasic(cbind(reference_main$meta_data, reference_main$umap$embedding), 
#          title = 'Reference Cells in Reference UMAP Space of Main celltypes', 
#          color.by = subtype_col_name,
#          save_path = paste(output_dir, '/1-02-sub_celltype_refer.png', sep=''),
#          size = c(8,10.5))



# Map main celltype
print('Mapping main celltypes...')
query_exp = seurat_merged$RNA$data
query_metadata= seurat_merged@meta.data
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference_main,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP
query = knnPredict(query, 
                   reference_main, 
                   reference_main$meta_data[,maintype_col_name], 
                   k = k,
                   save_as = paste(maintype_col_name,'.pred',sep=''))
plotBasic(cbind(query$meta_data, query$umap), 
          title = 'Query Cells in Reference UMAP Space', 
          color.by = paste(maintype_col_name,'.pred',sep=''),
          save_path = paste(output_dir, '/2-main_celltype_pred.png', sep=''),
          color.mapping = color.mapping.main,
          size = c(8,10))


########## Predict Sub Celltype ##########
print('Start working on sub celltypes.')
celltypes_with_sub = list()
celltypes_without_sub = list()
for (main_type in names(table(reference_main$meta_data[,maintype_col_name]))){
	ref_metadata_sub <- reference_main$meta_data %>% 
                        filter(get(maintype_col_name) == main_type) %>% 
                        #filter(get(subtype_col_name) != '') %>% 
                        rownames_to_column() %>% column_to_rownames()
	if (length(table(ref_metadata_sub[,subtype_col_name])) <= 1){
        celltypes_without_sub[main_type] <- names(table(ref_metadata_sub[,subtype_col_name]))[1]
    } else{
        celltypes_with_sub[main_type] <- as.data.frame(names(table(ref_metadata_sub[,subtype_col_name])))
    }
}
print("Main celltypes with sub celltypes:")
print(celltypes_with_sub)
print("Main celltypes without sub celltypes:")
print(celltypes_without_sub)

# plot distribution of celltypes in reference space
for (main_type in names(celltypes_with_sub)){
    ## load reference for main celltype
    reference_sub <- readRDS(paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds', sep=''))
    # If it cannot automatically find uwot model, manually set the path.
    #reference_sub$save_uwot_path <- ''

    # plot distribution of celltypes in reference space
    plotBasic(cbind(reference_sub$meta_data, reference_sub$umap$embedding), 
              title = paste('Reference Cells in Reference UMAP Space of', main_type), 
              color.by = subtype_col_name,
              color.mapping = color.mapping.sub,
              save_path = paste(output_dir, '/3-sub_celltype_refer_',gsub(' ','-',main_type),'.png', sep=''),
              size = c(8,10))
}

# Map sub celltype
print('Mapping sub celltypes...')
i = 1
sub_results = data.frame()
for (main_type in names(celltypes_with_sub)){
    reference_sub <- readRDS(paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds', sep=''))
	query_metadata_sub = query$meta_data[query$meta_data[paste(maintype_col_name,'.pred',sep='')] == main_type,]
	query_exp_sub = as.array(query$exp)[,rownames(query_metadata_sub)]
	query_sub = mapQuery(query_exp_sub,             # query gene expression (genes x cells)
	                     query_metadata_sub,        # query metadata (cells x attributes)
	                     reference_sub,             # Symphony reference object
	                     do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
	                     do_umap = TRUE)        # project query cells into reference UMAP
	query_sub = knnPredict(query_sub, 
	                       reference_sub, 
	                       reference_sub$meta_data[,subtype_col_name], 
	                       k = k,
	                       save_as = paste(main_type,'.pred',sep=''))
	umap_combined_labels = cbind(query_sub$meta_data, query_sub$umap)
    plotBasic(umap_combined_labels, 
              title = paste('Query Cells in Reference UMAP Space of', main_type), 
              color.by = paste(main_type,'.pred',sep=''),
              color.mapping = color.mapping.sub,
              save_path = paste(output_dir, '/4-sub_celltype_pred_',gsub(' ','-',main_type),'.png', sep=''),
              size = c(8,10))
    query_sub$meta_data$celltype.pred.combined <- query_sub$meta_data[,paste(main_type,'.pred',sep='')]
    query_sub$meta_data[paste(main_type,'.UMAP1',sep='')] <- query_sub$umap['UMAP1']
    query_sub$meta_data[paste(main_type,'.UMAP2',sep='')] <- query_sub$umap['UMAP2']
    if (i==1){
        sub_results = query_sub$meta_data
        i <- i+1
    } else {
        sub_results = bind_rows(sub_results, query_sub$meta_data)
    }
}

## merge results
# add umap location in main reference
query$meta_data <- cbind(query$meta_data, query$umap)
# process celltypes without subtype
query$meta_data$celltype.pred.combined <- query$meta_data[,paste(maintype_col_name,'.pred',sep='')]
for (main_type in names(celltypes_without_sub)){
    query$meta_data$celltype.pred.combined <- replace(query$meta_data$celltype.pred.combined, query$meta_data$celltype.pred.combined==main_type, celltypes_without_sub[main_type])
}
label_main <- query$meta_data %>% filter(!(celltype.pred.combined %in% names(celltypes_with_sub)))
print(str(label_main))
label_final <- bind_rows(label_main, sub_results)
label_final$celltype.pred.combined <- as.factor(as.character(label_final$celltype.pred.combined))

saveRDS(label_final, file=paste(output_dir, '/symphony_celltype_results.rds', sep=''))
write.csv(label_final, file=paste(output_dir, '/symphony_celltype_results.csv', sep=''), row.names = TRUE)

print(str(label_final))
print(table(label_final$celltype.pred.combined))


#############################################################################################
########## Draw Celltype Proportion ##########
ref_metadata_cleaned <- ref_metadata_for_sub #%>% filter(get(subtype_col_name) != '')
label_final <- readRDS(paste(output_dir, '/symphony_celltype_results.rds', sep=''))

for (p in c(output_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}
## main celltypes
# frequences for reference and queries
freq_refer <- as.data.frame(table(ref_metadata_cleaned[,maintype_col_name]))
freq_query <- as.data.frame(table(label_final[,paste(maintype_col_name,'.pred',sep='')]))
# add cell types with 0 counts
freq_refer$Var1 <- as.character(freq_refer$Var1)
freq_query$Var1 <- as.character(freq_query$Var1)
for (celltype in freq_refer$Var1){
    if (celltype %in% freq_query$Var1){next}
    freq_query = rbind(freq_query, data.frame(Var1=c(celltype), Freq=c(0)))
}
freq_refer = freq_refer %>% arrange(Var1)
freq_query = freq_query %>% arrange(Var1)

proportions.main <- data.frame(col1=names(table(ref_metadata_cleaned[,maintype_col_name])),
                               col2=freq_refer$Freq,
                               col3=freq_query$Freq,
                               col4=freq_refer$Freq/sum(freq_refer$Freq),
                               col5=freq_query$Freq/sum(freq_query$Freq),
                               check.rows=TRUE)
names(proportions.main) <- c(paste(maintype_col_name),
                             paste(maintype_col_name,'.freq.refer',sep=''),
                             paste(maintype_col_name,'.freq.pred',sep=''),
                             paste(maintype_col_name,'.prop.refer',sep=''),
                             paste(maintype_col_name,'.prop.pred',sep=''))
# save csv
write.csv(proportions.main, file=paste(output_dir, '/5-celltype_proportion_main.csv', sep=''), row.names = FALSE)

plotProp <- function(proportions,
                     celltype_col_name,     # the column you want to used as celltypes
                     x_col_name,            # the column for x
                     y_col_name,            # the column for y
                     x_label = 'Symphony Results',
                     y_label = 'Reference',
                     title = 'Proportions of Celltypes',         # Plot title
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right', # Show cell type legend
					 save_path = './plot.png',
                     size = c(8,10)) {  
    # set limits of axis
    axis_limit <- ceiling(max(max(proportions[,x_col_name]),max(proportions[,y_col_name])) * 10) / 10
    # calculate spearman corr
    spearman <- cor.test(proportions[,x_col_name], proportions[,y_col_name], method='spearman')
    
    # plot
    p = ggplot(proportions, aes(x=get(x_col_name), y=get(y_col_name))) + 
               geom_abline() + 
               geom_point_rast(aes(col = get(celltype_col_name)))
    if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    p = p + theme_bw() +
        labs(title = title, color = celltype_col_name) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(legend.position='right') +
        theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    	theme(aspect.ratio=1) +
        xlim(0, axis_limit) +
    	ylim(0, axis_limit) +
    	labs(x = x_label, y = y_label) +
        guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
    # add results of Spearman Corr.
    p = cowplot::add_sub(p, paste('Spearman Corr. =', spearman$estimate))
    p = cowplot::add_sub(p, paste('p =', spearman$p.value))

    # save plot
    ggsave(save_path, plot=p, width=size[2], height=size[1])

    return(spearman)
}


# plot main celltypes
spearman <- plotProp(proportions = proportions.main,
                     celltype_col_name = paste(maintype_col_name),
                     x_col_name = paste(maintype_col_name,'.prop.pred',sep=''),
                     y_col_name = paste(maintype_col_name,'.prop.refer',sep=''),
                     color.mapping = color.mapping.main,
                     save_path = paste(output_dir, '/5-celltype_proportion_main.png', sep=''))
# record spearman correlation
spearman.results <- data.frame(Item = c('Main Cell Types'),
                               Spearman.Correlation = c(spearman$estimate),
                               p.value = c(spearman$p.value))


## sub celltypes
freq_refer <- as.data.frame(table(ref_metadata_cleaned[,subtype_col_name]))
freq_query <- as.data.frame(table(label_final$celltype.pred.combined))
# add cell types with 0 counts
freq_refer$Var1 <- as.character(freq_refer$Var1)
freq_query$Var1 <- as.character(freq_query$Var1)
for (celltype in freq_refer$Var1){
    if (celltype %in% freq_query$Var1){next}
    freq_query = rbind(freq_query, data.frame(Var1=c(celltype), Freq=c(0)))
}
freq_refer = freq_refer %>% arrange(Var1)
freq_query = freq_query %>% arrange(Var1)

proportions.sub <- data.frame(col1=freq_refer$Var1,
                              col2=freq_refer$Freq,
                              col3=freq_query$Freq,
                              col4=freq_refer$Freq/sum(freq_refer$Freq),
                              col5=freq_query$Freq/sum(freq_query$Freq),
                              check.rows=TRUE)
names(proportions.sub) <- c(paste(subtype_col_name),
                            paste(subtype_col_name,'.freq.refer',sep=''),
                            paste(subtype_col_name,'.freq.pred',sep=''),
                            paste(subtype_col_name,'.prop.refer',sep=''),
                            paste(subtype_col_name,'.prop.pred',sep=''))
# save csv
write.csv(proportions.sub, file=paste(output_dir, '/5-celltype_proportion_sub.csv', sep=''), row.names = FALSE)
# plot all subtypes and save
spearman <- plotProp(proportions = proportions.sub,
                     celltype_col_name = paste(subtype_col_name),
                     x_col_name = paste(subtype_col_name,'.prop.pred',sep=''),
                     y_col_name = paste(subtype_col_name,'.prop.refer',sep=''),
                     color.mapping = color.mapping.sub,
                     save_path = paste(output_dir, '/5-celltype_proportion_sub_all.png', sep=''))
# record spearman correlation
spearman.results <- rbind(spearman.results, list('All Sub-celltypes', spearman$estimate, spearman$p.value))


# plot by main celltypes
for (main_type in names(table(ref_metadata_cleaned[,maintype_col_name]))){
    # choose cells only in this main celltype form the metadata
    ref_metadata_sub <- ref_metadata_cleaned %>% filter(get(maintype_col_name) == main_type) #%>% filter(get(subtype_col_name) != '') 
    # check number of sub celltype, if <= 1, skip it.
    if (length(table(ref_metadata_sub[,subtype_col_name])) <= 1){next}
	
    # choose proportions of subtypes within the main celltypes
    subgroup = proportions.sub %>% filter(get(subtype_col_name) %in% names(table(ref_metadata_sub[,subtype_col_name])))
    # plot and calculate spearman correlation
    spearman <-plotProp(proportions = subgroup,
                        celltype_col_name = paste(subtype_col_name),
                        x_col_name = paste(subtype_col_name,'.prop.pred',sep=''),
                        y_col_name = paste(subtype_col_name,'.prop.refer',sep=''),
                        title = paste('Proportions of Subtypes of', main_type), 
                        color.mapping = color.mapping.sub,
                        save_path = paste(output_dir, '/5-celltype_proportion_sub_',gsub(' ','-',main_type),'.png', sep=''))
    spearman.results <- rbind(spearman.results, list(main_type, spearman$estimate, spearman$p.value))
}

write.csv(spearman.results, file=paste(output_dir, '/5-celltype_proportion_spearman.csv', sep=''), row.names = FALSE)


########## Draw Bubble Plot ##########
count <- as.matrix(seurat_merged$RNA$counts)
label_final <- readRDS(paste(output_dir, '/symphony_celltype_results.rds', sep=''))

subtype <- data.frame(cellid = rownames(label_final), celltype=label_final[,paste(maintype_col_name,'.pred',sep='')])
subtype$celltype <- as.character(subtype$celltype)
count <- count[,subtype$cellid]

metadata <- data.frame(cell_id = subtype$cellid, cell_types = subtype$celltype)
rownames(metadata) <- metadata$`cell_id`
seurat <- CreateSeuratObject(counts = count, meta.data = metadata)
seurat <- NormalizeData(seurat)

plotBubble <- function(seurat_obj,
                       title = 'Bubble Plot',         # Plot title
					   save_path = './plot.png',
                       size = c(8,10)) { 
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat@meta.data$cell_types <- as.factor(seurat@meta.data$cell_types)
Idents(seurat) <- seurat@meta.data$cell_types

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extracting top 5 marker genes for each cell type, sorted by avg_log2FC in descending order
top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Creating a new column to rank genes within each cell type, sorting from high to low
top5_markers <- top5_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  mutate(gene_rank = row_number())

# Formatting the data for the bubble plot
# Extracting log fold change and proportion of cells expressing the marker genes
bubble_data <- top5_markers %>% 
  select(gene, cell_type = cluster, logFC = avg_log2FC, proportion = pct.1, gene_rank)

# Sorting genes for left-diagonal layout
bubble_data <- bubble_data %>%
  arrange(desc(cell_type), desc(gene_rank))

# Converting gene names to factors to maintain the order in the plot
bubble_data$gene <- factor(bubble_data$gene, levels = rev(unique(bubble_data$gene)))

# Creating the bubble plot
bubble_plot <- ggplot(bubble_data, aes(x = cell_type, y = gene, size = proportion, color = logFC)) +
  geom_point() +
  scale_color_gradient2(low = "white", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = title,
       x = "Cell Type",
       y = "Gene",
       size = "Proportion of Cells Expressing",
       color = "Log Fold Change")

# Save the bubble plot
ggsave(save_path, plot = bubble_plot, width = size[2], height = size[1])
                     }

# for main celltypes
plotBubble(seurat_obj=seurat,
           title="Marker Genes of Main Celltypes",
           save_path=paste(output_dir, '/6-bubble_plot_main.png', sep=''))

## for sub-types
count <- as.matrix(seurat_merged$RNA$data)
label_final <- readRDS(paste(output_dir, '/symphony_celltype_results.rds', sep=''))

subtype <- data.frame(cellid = rownames(label_final), celltype=label_final$celltype.pred.combined)
subtype$celltype <- as.character(subtype$celltype)
count <- count[,subtype$cellid]
# plot by main celltypes
for (main_type in names(table(label_final[,paste(maintype_col_name,'.pred',sep='')]))){
    # choose cells only in this main celltype form the metadata
    label_final_sub <- label_final %>% 
                       filter(get(paste(maintype_col_name,'.pred',sep='')) == main_type) #%>% 
                       #filter(celltype.pred.combined != '') 
    label_final_sub$celltype.pred.combined <- as.character(label_final_sub$celltype.pred.combined)
    # check number of sub celltype, if <= 1, skip it.
    if (length(table(label_final_sub$celltype.pred.combined)) <= 1){next}

    # choose proportions of subtypes within the main celltypes
    subtype_sub <- subtype[subtype$celltype %in% names(table(label_final_sub$celltype.pred.combined)),]
    count_sub <- count[,subtype_sub$cellid]

    metadata <- data.frame(cell_id = subtype_sub$cellid, cell_types = subtype_sub$celltype)
    rownames(metadata) <- metadata$`cell_id`
    seurat <- CreateSeuratObject(counts = count_sub, meta.data = metadata)
    seurat <- NormalizeData(seurat)

    # plot 
    plotBubble(seurat_obj=seurat,
               title=paste("Marker Genes of", main_type),
               save_path=paste(output_dir, '/6-bubble_plot_',gsub(' ','-',main_type),'.png', sep=''))
}
