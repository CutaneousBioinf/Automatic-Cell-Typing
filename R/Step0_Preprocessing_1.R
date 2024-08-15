# This script is used to preprocess the Xenium data and generate the Giotto object

args = commandArgs(trailingOnly=TRUE) 

if (length(args)<2){
    print("Usage: xenium.folders results.folder")
    quit(status = 1)
}

xenium_folders <- scan(args[1],what="",quiet=T)
results_folder <- args[2]

# results_folder = "/home/zhanh/Janssen_AD_xenium/results"
# xenium_folders=scan(paste0("Janssen_AD_xenium/Janssen_Xenium_AD_UMcustom.panel"),what="",quiet=T)


#### check and set up Giotto environment ####
### load the latest Matrix
#library(Matrix, lib.loc = "/usr/local/lib/R/moreLibs")
library(Matrix)

# Ensure Giotto Suite is installed.
if (!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

library(Giotto)
# check if all following packages are installed, if not then install
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,RColorBrewer,ggplot2,dplyr,hrbrthemes,viridis,ggpubr,sjmisc) 

# Ensure the Python environment for Giotto has been installed.
genv_exists <- checkGiottoEnvironment()
if (!genv_exists) {
  #The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

#color generator: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# it generates 74 colors for plots
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# set working directory
#results_folder = '/path/to/save/directory/'
#results_folder = getwd()
setwd(results_folder)

# set giotto python path
# set python path to your preferred python version path
# set python path to NULL if you want to use the giotto miniconda environment
python_path = NULL

# Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)

xenium_gobj_list=list()
raw_tx_profile=data.frame()

for(xfolder in xenium_folders){ ##MTP: Since the number of folders varies by tissue, I use a loop - this also makes the code tidier
    # general files (some are supplemental files)
    #settings_path = paste0(xfolder, 'experiment.xenium') ##MTP: commented out, as does not seem to be used

    # files (SUBCELLULAR):
    cell_bound_path = paste0(xfolder, '/cell_boundaries.csv.gz')
    nuc_bound_path = paste0(xfolder, '/nucleus_boundaries.csv.gz')
    tx_path = paste0(xfolder, '/transcripts.csv.gz')
    feat_meta_path = paste0(xfolder, '/cell_feature_matrix/features.tsv.gz') # (also used in aggregate)

    # files (AGGREGATE):
    #expr_mat_path = paste0(xfolder, 'cell_feature_matrix') ##MTP: commented out, as does not seem to be used
    #cell_meta_path = paste0(xfolder, 'cells.csv.gz') # contains spatlocs ##MTP: commented out, as does not seem to be used

    #### Xenium feature types exploration ####
    # load features metadata
    # (make sure cell_feature_matrix folder is unpacked)
    feature_dt = data.table::fread(feat_meta_path, header = FALSE)
    colnames(feature_dt) = c('feat_ID','feat_name','feat_type')

    # find the feature IDs that belong to each feature type
    feat_types = names(feature_dt[, table(feat_type)])
    feat_types_IDs = lapply(
    feat_types, function(type) feature_dt[feat_type == type, unique(feat_name)] # use feat_name instead of feat_ID here
    )
    names(feat_types_IDs) = feat_types

    #### loading Xenium data and create Giotto object ####
    # load transcript-level data
    tx_dt = data.table::fread(tx_path)
    data.table::setnames(x = tx_dt,
                        old = c('feature_name', 'x_location', 'y_location'),
                        new = c('feat_ID', 'x', 'y'))
    
    
    ### AT: NOTE Cosmx/Xenium definition of unassigned transcripts are different
    tx_dt = tx_dt[cell_id > 0 & cell_id!="UNASSIGNED"]
    cat('Transcripts info available:\n ', paste0('"', colnames(tx_dt), '"'), '\n',
        'with', tx_dt[,.N], 'unfiltered detections\n')
    raw_tx_profile=rbind(raw_tx_profile,tx_dt) ##MTP: Creating this in the loop, rather than at the end, for simplicity

    # filter by qv (Phred score)
    tx_dt_filtered = tx_dt[qv >= 20]
    cat('and', tx_dt_filtered[,.N], 'filtered detections\n\n')

    # separate detections by feature type
    tx_dt_types = lapply(
    feat_types_IDs, function(types) tx_dt_filtered[feat_ID %in% types]
    )

    invisible(lapply(seq_along(tx_dt_types), function(x) {
    cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][,.N], '\n')
    }))

    gpoints_list = lapply(
    tx_dt_types, function(x) createGiottoPoints(x = x)
    )
    
    # Load polygon data
    cellPoly_dt = data.table::fread(cell_bound_path)
    nucPoly_dt = data.table::fread(nuc_bound_path)

    data.table::setnames(cellPoly_dt,
                        old = c('cell_id', 'vertex_x', 'vertex_y'),
                        new = c('poly_ID', 'x', 'y'))

    data.table::setnames(nucPoly_dt,
                        old = c('cell_id', 'vertex_x', 'vertex_y'),
                        new = c('poly_ID', 'x', 'y'))

    # only keep cells with valid assigned transcripts
    cellPoly_dt = cellPoly_dt[poly_ID %in% unique(tx_dt_filtered$cell_id)]
    nucPoly_dt = nucPoly_dt[poly_ID %in% unique(tx_dt_filtered$cell_id)]

    gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                            name = 'cell',
                                            calc_centroids = TRUE)

    gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                            name = 'nucleus',
                                            calc_centroids = TRUE)

    # Create Giotto Object
    gobj = createGiottoObjectSubcellular(
    gpoints = list(rna = gpoints_list$`Gene Expression`,
                    ### blank_code = gpoints_list$`Unassigned Codeword`,
                    neg_code = gpoints_list$`Negative Control Codeword`,
                    neg_probe = gpoints_list$`Negative Control Probe`),
    gpolygons = list(cell = gpoly_cells,
                    nucleus = gpoly_nucs),
    instructions = instrs
    )

    ## generate aggregated expression based on feature and boundary (polygon) information
    # Calculate the overlaps of the 'rna' feature data within the 'cell' polygon boundary info.
    # Find the feature points overlapped by polygons. This overlap information is then
    # returned to the relevant giottoPolygon object's overlaps slot

    gobj = calculateOverlapRaster(gobj,spatial_info = 'cell', feat_info = 'rna')

    gobj = calculateOverlapRaster(gobj, spatial_info = 'cell', feat_info = 'neg_probe')
    gobj = calculateOverlapRaster(gobj, spatial_info = 'cell', feat_info = 'neg_code')

    ### can check through get_polygon_info(fov_join, polygon_name="cell",polygon_overlap="rna")

    # Convert the overlap information into a cell by feature expression matrix which
    # is then stored in the Giotto object's expression slot


    # Assign polygon overlaps information to expression matrix
    gobj = overlapToMatrix(gobj,feat_info = 'rna',poly_info = 'cell',name = 'raw')
    gobj = overlapToMatrix(gobj, feat_info = 'neg_probe',poly_info = 'cell',name = 'raw')
    gobj = overlapToMatrix(gobj, feat_info = 'neg_code',poly_info = 'cell',name = 'raw')


    # get FOV info from loaded transcript-level data
    fov_info = as.data.frame(tx_dt_filtered[,c("cell_id","fov_name")])
    fov_info=fov_info[fov_info$cell_id!="UNASSIGNED",] ##MTP: Somehow there are some unsassigned transcripts 
    fov_info = fov_info[!duplicated(fov_info$cell_id),]
    fov_info = fov_info[order(fov_info$cell_id, decreasing = F),]

    # add FOV info to cell metadata
    gobj = addCellMetadata(gobj, new_metadata = fov_info$fov_name, vector_name = "fov")

    xenium_gobj_list=append(xenium_gobj_list,gobj) ##MTP: Create list of gobjs to be joined
}





############################################################################
# if you need to combine giotto objects loaded from different CosMx directories into a single object
# below is an example of two objects, will extend to >2 objects
obj_ids = paste0("obj",1:length(xenium_gobj_list))
if(length(xenium_gobj_list)>1){ ##MTP: Important not to use join when there is only one gobj
xenium_gobj = joinGiottoObjects(gobject_list = xenium_gobj_list,
                                gobject_names = obj_ids,
                                x_padding = 8000, y_padding=2000)
} else { xenium_gobj=xenium_gobj_list[[1]] }



#### visualize Giotto object and cells ####
spatInSituPlotPoints(xenium_gobj,show_image = F,image_name = NULL,spat_unit = 'cell',point_size = 0.01,show_polygon = TRUE,use_overlap = FALSE,polygon_feat_type = 'cell',polygon_color = 'white',polygon_line_size = 0.03,coord_fix_ratio = TRUE,save_param = list(base_height = 10, dpi = 1000,save_name = '1_inSituFovs'))

# Plot the generated centroids information
spatPlot2D(xenium_gobj,
        spat_unit = 'cell',
        point_shape = 'no_border',
        point_size = 0.1,
        point_alpha = 0.4,
        save_param = list(
            base_width = 10,
            base_height = 10,
            save_name = '2_spatplot'))

###################################################
### QC plots
calulate_per_fov <- function(slot){
if (slot == "raw"){
    counts = as.matrix(xenium_gobj@expression$cell$rna$raw@exprMat)
} else {
    counts = as.matrix(xenium_gobj@expression$cell$rna$normalized@exprMat)
}

meta = pDataDT(xenium_gobj)
out = sapply(unique(meta$fov), function(x){
    cellIDs = subset(meta,fov==x)$cell_ID
    raw_tx_profile = as.data.frame(raw_tx_profile)
    n_tx = nrow(raw_tx_profile[raw_tx_profile$fov_name==x,])
    data.frame(fov=x,
            tx_assigned = sum(counts[,cellIDs])/n_tx,
            count_tx = sum(counts[,cellIDs]),
            count_cells =length(cellIDs),
            tx_cell =sum(counts[,cellIDs])/length(cellIDs))
}, simplify=F) %>% bind_rows()
return(out)
}

### raw
## At FOV level
tb = calulate_per_fov("raw")

colnames(tb) = c("fov","Transcript assigned rate","Total number of transcripts","Total number of cells","Transcript-cell ratio")

pdf("2.1.a_QC_FOV_level.pdf",heigh=10,width=10)
###  width = 6000, height = 6000, res = 600
    p1 = ggplot(tb, aes(x=fov, y=`Transcript assigned rate`)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(limits = c(0,1))
    p2 = ggplot(tb, aes(x=fov, y=`Total number of transcripts`)) + geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p3 = ggplot(tb, aes(x=fov, y=`Total number of cells`)) + geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p4 = ggplot(tb, aes(x=fov, y=`Transcript-cell ratio`)) + geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p5 = ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)
    annotate_figure(p5, top = text_grob("QC metrics at FOV level", face = "bold", size = 14))
dev.off()

    write.table(tb, "2.1.b_QC_FOV_level.csv", row.names = F, col.names = T, sep=",")

## At cell level
tb_cell = pDataDT(xenium_gobj)
counts = as.matrix(xenium_gobj@expression$cell$rna$raw@exprMat)
tb_cell$count_tx = colSums(counts)
tb_cell$count_gene = colSums(counts!=0)
tb_cell = as.data.frame(tb_cell)
colnames(tb_cell)[colnames(tb_cell)=="nr_feats" | colnames(tb_cell)=="count_gene"]  <- "Number of detected genes per cell"
colnames(tb_cell)[colnames(tb_cell)=="total_expr" | colnames(tb_cell)=="count_tx"]  <- "Number of transcripts per cell"

pdf("2.2.a_QC_cell_level.pdf", width=10,height=10)
#width = 6000, height = 3000, res = 600)
p1 = ggplot(tb_cell, aes(x=fov, y=`Number of transcripts per cell`)) + geom_boxplot() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 = ggplot(tb_cell, aes(x=fov, y=`Number of detected genes per cell`)) + geom_boxplot() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 = ggarrange(p1,p2,ncol = 2,nrow = 1)
annotate_figure(p3, top = text_grob("QC metrics at cell level", face = "bold", size = 14))
dev.off()

write.table(tb_cell, "2.2.b_QC_cell_level.csv", row.names = F, col.names = T, sep=",")

###################################
## change of cell counts under different gene/transcript cutoffs
pre_tb = function(df,col){
    df = df[order(df[,col]),]
    tb = table(df[,col])
    tb = as.data.frame(cbind(cutoff=as.numeric(names(tb)), count=cumsum(tb)))
    tb$pc = 1-tb$count/nrow(df)
    df_out = tb[,-2]
    df_out$fov = "all"
    for (fov in sort(unique(df$fov))){
        df_sub = df[df$fov==fov,]
        tb = table(df_sub[,col])
        tb = as.data.frame(cbind(cutoff=as.numeric(names(tb)), count=cumsum(tb)))
        tb$pc = 1-tb$count/nrow(df_sub)
        tb = tb[,-2]
        tb$fov = fov
        df_out = rbind(df_out,tb)
    }
    return(df_out)
}
# cutoff of transcripts
tb_cutoff1 = pre_tb(tb_cell,"Number of transcripts per cell")
tb_cutoff2 = pre_tb(tb_cell,"Number of detected genes per cell")
tx_cutoff = max(tb_cutoff1[tb_cutoff1$fov=="all"&tb_cutoff1$pc>0.95,"cutoff"])
gene_cutoff = max(tb_cutoff2[tb_cutoff2$fov=="all"&tb_cutoff2$pc>0.95,"cutoff"])


pdf("2.3_cell_count_change.pdf",width=20,height=10)
###, width = 10000, height = 4000, res = 500)
p1 = ggplot(tb_cutoff1, aes(x=cutoff, y=pc, color=fov)) + theme(text = element_text(size = 10)) + 
scale_x_continuous(breaks = sort(c(seq(0,max(tb_cutoff1$cutoff),50),tx_cutoff))) + 
scale_y_continuous(breaks = sort(c(seq(0,1,0.2),0.95))) + geom_point() + 
xlab("Threshold of transcripts #") + ylab("Proportion of cells retained") + 
geom_vline(xintercept=tx_cutoff,linetype=2) + geom_hline(yintercept=0.95,linetype=2)
p2 = ggplot(tb_cutoff2, aes(x=cutoff, y=pc, color=fov)) + theme(text = element_text(size = 10)) + scale_x_continuous(breaks = sort(c(seq(0,max(tb_cutoff2$cutoff),50),gene_cutoff))) + scale_y_continuous(breaks = sort(c(seq(0,1,0.2),0.95))) + geom_point() + 
xlab("Threshold of detected genes #") + ylab("Proportion of cells retained") + 
geom_vline(xintercept=gene_cutoff,linetype=2) + geom_hline(yintercept=0.95,linetype=2)
p3 = ggarrange(p1,p2,ncol = 1,nrow = 2,common.legend = TRUE, legend="right")
annotate_figure(p3, top = text_grob("Change of cell counts under different thresholds of transcripts/detected genes #", face = "bold", size = 14), bottom = text_grob("*The vertical dashed line indicates the maximum threshold to retain >95% of cells across all FOVs", size = 14))
dev.off()

############################################
# Data filtering
xenium_gobj = filterGiotto(gobject = xenium_gobj,
                        spat_unit = 'cell',
                        poly_info = 'cell',
                        expression_threshold = 1,
                        feat_det_in_min_cells = 5,
                        min_det_feats_per_cell = gene_cutoff)

# Add data statistics
xenium_gobj = addStatistics(xenium_gobj, expression_values = 'raw')

#######################################################
### after filtering

## QC plots after filtering
tb2 = calulate_per_fov("raw")
colnames(tb2) = c("fov","Transcript assigned rate","Total number of transcripts","Total number of cells","Transcript-cell ratio")
tb$Filter = "NA"
tb2$Filter = "Applied"
tb = rbind(tb,tb2)
tb$Filter = factor(tb$Filter, levels = c("NA", "Applied"))

# generate histograms
pdf("3.1.a_QC_FOV_level_filtered.pdf",width=10,height=10)
### width = 6000, height = 6000, res = 600)
p1 = ggplot(tb, aes(x=fov, y=`Transcript assigned rate`, fill=Filter)) + geom_bar(position="dodge", stat = "identity") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(limits = c(0,1))
p2 = ggplot(tb, aes(x=fov, y=`Total number of transcripts`, fill=Filter)) + geom_bar(position="dodge", stat = "identity") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 = ggplot(tb, aes(x=fov, y=`Total number of cells`, fill=Filter)) + geom_bar(position="dodge", stat = "identity") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 = ggplot(tb, aes(x=fov, y=`Transcript-cell ratio`, fill=Filter)) + geom_bar(position="dodge", stat = "identity") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p5 = ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = TRUE,legend="bottom")
annotate_figure(p5, top = text_grob("QC metrics at FOV level after filtering", face = "bold", size = 14))
dev.off()

write.table(tb2, "3.1.b_QC_FOV_level_filtered.csv", row.names = F, col.names = T, sep=",")

## At cell level
tb_cell2 = pDataDT(xenium_gobj)
counts = as.matrix(xenium_gobj@expression$cell$rna$raw@exprMat)
tb_cell2$count_tx = colSums(counts)
tb_cell2$count_gene = colSums(counts!=0)
tb_cell2 = as.data.frame(tb_cell2)
colnames(tb_cell2)[colnames(tb_cell2)=="nr_feats"]  <- "Number of detected genes per cell"
colnames(tb_cell2)[colnames(tb_cell2)=="total_expr"]  <- "Number of transcripts per cell"

tb_cell$Filter = "NA"
tb_cell2$Filter = "Applied"
tb_cell = rbind(tb_cell[,colnames(tb_cell)],tb_cell2[,colnames(tb_cell)])
tb_cell$Filter = factor(tb_cell$Filter, levels = c("NA", "Applied"))


pdf("3.2.a_QC_cell_level_filtered.pdf",height=5,width=10)

### width = 6000, height = 3000, res = 600)
p1 = ggplot(tb_cell, aes(x=fov, y=`Number of transcripts per cell`, fill=Filter)) + geom_boxplot() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 = ggplot(tb_cell, aes(x=fov, y=`Number of detected genes per cell`, fill=Filter)) + geom_boxplot() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 = ggarrange(p1,p2,ncol = 2,nrow = 1,common.legend = TRUE,legend="bottom")
annotate_figure(p3, top = text_grob("QC metrics at cell level after filtering", face = "bold", size = 14))
dev.off()

write.table(tb_cell2, "3.2.b_QC_cell_level_filtered.csv", row.names = F, col.names = T, sep=",")

###############################

# Normalize expression
xenium_gobj = normalizeGiotto(gobject = xenium_gobj,
                            feat_type='rna',norm_methods='standard',spat_unit = 'cell',
                            verbose = T)


# standard method of normalization (log normalization based)
xenium_gobj <- normalizeGiotto(gobject = xenium_gobj,
                            feat_type = 'neg_probe',
                            norm_methods = 'standard',
                            spat_unit = 'cell',
                            library_size_norm = FALSE,
                            verbose = TRUE)

xenium_gobj <- normalizeGiotto(gobject = xenium_gobj,
                            feat_type = 'neg_code',
                            norm_methods = 'standard',
                            spat_unit = 'cell',
                            library_size_norm = FALSE,
                            verbose = TRUE)

# add statistics based on log normalized values for features rna and negative probes
xenium_gobj = addStatistics(gobject = xenium_gobj, expression_values = 'normalized',feat_type = 'rna')


# QC plots after normalization
tb3 = calulate_per_fov("normalized")

colnames(tb3) = c("fov","Transcript assigned rate","Total number of normalized transcripts","Total number of cells","Normalized transcript-cell ratio")

tb_cell3 = pDataDT(xenium_gobj)
colnames(tb_cell3)[c(2,5)] = c("fov","Normalized transcript count per cell")

# generate plots
pdf("4.1_QC_FOV_level_normalized.pdf",width=15,height=5)
### width = 9000, height = 3000, res = 600)
p1 = ggplot(tb_cell3, aes(x=`Normalized transcript count per cell`)) + geom_histogram(binwidth=25)
p2 = ggplot(tb_cell3, aes(x=fov, y=`Normalized transcript count per cell`)) + geom_boxplot() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 = ggplot(tb3, aes(x=fov, y=`Normalized transcript-cell ratio`)) + geom_bar(stat = "identity") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 = ggarrange(p1,p2,p3,ncol = 3,nrow = 1)
annotate_figure(p4, top = text_grob("QC metrics after normalization", face = "bold", size = 14))
dev.off()

####################################################
# Plot spatially as centroids
spatPlot2D(gobject = xenium_gobj,
        cell_color = 'total_expr',
        color_as_factor = FALSE,
        show_image = F,
        image_name = NULL,
        point_size = 0.25,
        point_alpha = 0.75,
        save_param = list(base_height = 10, dpi = 1000,
                            save_name = '4.2_NormExp_centroids'))


# Plot spatially as color-scaled polygons
spatInSituPlotPoints(xenium_gobj,
                    show_polygon = TRUE,
                    polygon_color = 'gray',
                    polygon_line_size = 0.03,
                    polygon_fill = 'total_expr',
                    polygon_fill_as_factor = FALSE,
                    save_param = list(base_height = 10, dpi = 1000,
                                    save_name = '4.3_NormExp_rna_polys'))


###################################################
# Calculate highly variable features
xenium_gobj = calculateHVF(gobject = xenium_gobj,
                        spat_unit = 'cell',
                        save_param = list(save_name = '5.1_HVF'))

cat(fDataDT(xenium_gobj)[, sum(hvf == 'yes')], 'hvf found', '\n')

#### dimension reduction and clustering ####
# PCA
xenium_gobj = runPCA(gobject = xenium_gobj,spat_unit = 'cell',expression_values = 'scaled',feats_to_use = NULL,scale_unit = F,center = F)


# Visualize Screeplot and PCA
screePlot(xenium_gobj,
        ncp = 20,
        save_param = list(save_name = '5.2_screePlot'))

plotPCA(xenium_gobj,
        spat_unit = 'cell',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2,
        save_param = list(save_name = '5.3_PCA'))

# UMAP
xenium_gobj = runUMAP(xenium_gobj,
                      dimensions_to_use = 1:10,
                      spat_unit = 'cell', dim_reduction_to_use="pca")

plotUMAP(xenium_gobj,
         point_size = 0.01,
         save_param = list(save_name = '5.4_UMAP'))

pc_num <-25

# run Harmony
xenium_gobj = runGiottoHarmony(xenium_gobj, 
                            vars_use = "fov",
                            dim_reduction_to_use = "pca",
                            dimensions_to_use = 1:pc_num)

# Generate UMAP from Harmony
xenium_gobj = runUMAP(xenium_gobj,
                    dim_reduction_to_use = "harmony",
                    dim_reduction_name = "harmony",
                    dimensions_to_use = 1:pc_num,
                    n_threads = 4)

plotUMAP(gobject = xenium_gobj, save_param = list(save_name = '5.5_UMAP_after_Harmony'))

# sNN and Leiden clustering
xenium_gobj = createNearestNetwork(xenium_gobj,
                                type = "sNN",
                                dim_reduction_to_use = "umap",
                                dim_reduction_name = "umap",
                                dimensions_to_use = 1:2, name="sNN.umap")

# the resolution here defaultly to be 0.1, but can be customized by the user
# the higher the resolution is, the more clusters there will be
xenium_gobj = doLeidenCluster(xenium_gobj,
                            network_name = "sNN.umap",
                            resolution = 0.05,
                            n_iterations = 1000)


# visualize UMAP cluster results
plotUMAP(gobject = xenium_gobj,
        cell_color = 'leiden_clus',
        cell_color_code = col_vector,
        show_NN_network = F,
        point_size = 1,
        save_param = list(save_name = '5.6_UMAP_leiden'))

# visualize UMAP cluster results, by batch
plotUMAP(gobject = xenium_gobj,
        cell_color = 'fov',
        cell_color_code = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
        show_NN_network = F,
        point_size = 1,
        save_param = list(save_name = '5.6_UMAP_byfov'))

# visualize UMAP cluster results, by batch
plotUMAP(gobject = xenium_gobj,
        cell_color = 'list_ID',
        cell_color_code = col_vector,
        ###cell_color_code = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
        show_NN_network = F,
        point_size = 1,
        save_param = list(save_name = '5.6_UMAP_bylist_ID'))
        
# Visualize UMAP and spatial results
spatInSituPlotPoints(xenium_gobj,
                    show_image = FALSE,
                    feats = NULL,
                    point_size = 0.05,
                    show_polygon = TRUE,
                    polygon_feat_type = 'cell',
                    polygon_alpha = 1,
                    polygon_color = 'black',
                    polygon_line_size = 0.01,
                    polygon_fill = 'leiden_clus',
                    polygon_fill_as_factor = TRUE,
                    polygon_fill_code = col_vector,
                    coord_fix_ratio = TRUE,
                    save_para = list(dpi = 1000,
                                    save_name = '5.6_insitu_leiden'))


spatInSituPlotPoints(xenium_gobj,
                    show_image = FALSE,
                    feats = NULL,
                    point_size = 0.05,
                    show_polygon = TRUE,
                    polygon_feat_type = 'cell',
                    polygon_alpha = 1,
                    polygon_color = 'black',
                    polygon_line_size = 0.01,
                    polygon_fill = 'fov',
                    polygon_fill_as_factor = TRUE,
                    polygon_fill_code = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
                    coord_fix_ratio = TRUE,
                    save_para = list(dpi = 1000,save_name = '5.6_insitu_FOV'))


spatInSituPlotPoints(xenium_gobj,
                    show_image = FALSE,
                    feats = NULL,
                    point_size = 0.05,
                    show_polygon = TRUE,
                    polygon_feat_type = 'cell',
                    polygon_alpha = 1,
                    polygon_color = 'black',
                    polygon_line_size = 0.01,
                    polygon_fill = 'list_ID',
                    polygon_fill_as_factor = TRUE,
                    polygon_fill_code = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
                    coord_fix_ratio = TRUE,
                    save_para = list(dpi = 1000,save_name = '5.6_insitu_list_ID'))

saveGiotto(xenium_gobj,"gobj",overwrite=T) ##MTP: We cannot just save the R state, because there is some data in C
