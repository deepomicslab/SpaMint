# library(SpaTalk)
options(warn = 0)  
# load starmap data
args = commandArgs(T)
st_dir <- args[1]
st_meta_dir <- args[2]
sc_coord_dir <- args[3]
meta_key <- args[4]
species <- args[5]
out_f <- args[6]
out_f = paste0(out_f,'/spa/')
print(length(args))
print(args)

if (length(args) > 6){
    n_cores = strtoi(args[7])
}else{
    n_cores = 4
}
# print(n_cores)

# TODO
args <- commandArgs(trailingOnly = FALSE)
scriptPath <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)]))
scriptPath <- dirname(scriptPath)
##########
dir.create(file.path(out_f), showWarnings = FALSE)
# print(out_f)
# print(meta_key)
if (grepl('csv', st_dir)){
    st_data = t(read.table(file = st_dir, sep = ',', header = TRUE,row.names = 1))
} else{
    st_data = t(read.table(file = st_dir, sep = '\t', header = TRUE,row.names = 1))
}

if (grepl('csv', st_meta_dir)){
    st_meta = read.table(file = st_meta_dir, sep = ',', header = TRUE,row.names = 1)
} else{
    st_meta = read.table(file = st_meta_dir, sep = '\t', header = TRUE,row.names = 1)
}

if (grepl('csv', sc_coord_dir)){
    sc_coord = read.table(file = sc_coord_dir, sep = ',', header = TRUE,row.names = 1)
} else{
    sc_coord = read.table(file = sc_coord_dir, sep = '\t', header = TRUE,row.names = 1)
}


if (species == 'Mouse'){
    if (grepl('spex', sc_coord_dir)){
        sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
    }else if (grepl('before', sc_coord_dir)){
        sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
    }else{
        sc_coord = sc_coord[c('x','y')]
    }
}

if (species == 'Human'){
    if (grepl('spex', sc_coord_dir)){
        sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
    }else{
        sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
    }
}

# if (grepl('mela', sc_coord_dir)){
#     if (grepl('spex', sc_coord_dir)){
#         sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
#     }else{
#         sc_coord = sc_coord[c('adj_spex_UMAP1','adj_spex_UMAP2')]
#     }
# }
print('loaded')

if (grepl('Human', species)){
    max_hop = 3
}else if ( 
    grepl('Mouse', species)){
    max_hop = 4
    # add neuronchat db
    lr_df_dir = paste0(scriptPath,"/../LR/mouse_LR_pairs.txt")
    new_lr_df = read.table(file = lr_df_dir, sep = '\t', header = FALSE)
    new_lr_df = new_lr_df[,c('V1','V2')]
    colnames(new_lr_df) = c('ligand','receptor')
    species = 'Mouse'
    new_lr_df$species = species
    new_lrpairs = rbind(lrpairs,new_lr_df)
    new_lrpairs$ligand = gsub("_", "-", new_lrpairs$ligand)
    new_lrpairs$receptor = gsub("_", "-", new_lrpairs$receptor)
    new_lrpairs = unique(new_lrpairs)
}



# subset by meta index
st_data = st_data[,rownames(sc_coord)]
# Formating
sc_coord$cell = rownames(sc_coord)
sc_coord$cell <- sub("^", "C",sc_coord$cell)
colnames(sc_coord) = c('x','y','cell')
sc_coord = sc_coord[,c('cell','x','y')]

colnames(st_data) = sc_coord$cell
colnames(st_data) = gsub("_", "-", colnames(st_data))
rownames(st_data) = gsub("_", "-", rownames(st_data))
st_data = as.data.frame(st_data)

obj <- createSpaTalk(st_data = as.matrix(st_data),
                     st_meta = sc_coord,
                     species = species,
                     if_st_is_sc = T,
                     spot_max_cell = 1,celltype = st_meta[[meta_key]])
tp_lst = unique(obj@meta$rawmeta$celltype)


obj <- find_lr_path(object = obj , lrpairs = new_lrpairs, pathways = pathways, if_doParallel = T, use_n_cores=n_cores, max_hop = max_hop)

for (tp1 in tp_lst) {
  for (tp2 in tp_lst) {
    if (tp1 != tp2) {
      tryCatch({
        obj <- dec_cci(object = obj, celltype_sender = tp1, celltype_receiver = tp2,
                       if_doParallel = T, use_n_cores = n_cores, pvalue = 0.1, n_neighbor = 20,
                       co_exp_ratio = 0.05, min_pairs = 2)
        obj <- dec_cci(object = obj, celltype_sender = tp2, celltype_receiver = tp1,
                       if_doParallel = T, use_n_cores = n_cores, pvalue = 0.1, n_neighbor = 20,
                       co_exp_ratio = 0.05, min_pairs = 2)
        print(tp1)
        print(tp2)
      }, error = function(e) {
        cat("Error occurred during iteration: tp1:", tp1, "tp2:", tp2, "Error:", conditionMessage(e), "\n")
      })
      next
    }
  }
}

# obj <- dec_cci_all(object = obj, if_doParallel = T, use_n_cores=n_cores, pvalue=0.1, n_neighbor = 20, co_exp_ratio=0.05,min_pairs=2)
write.csv(obj@lrpair, paste0(out_f,"/lr_pair.csv"), row.names = TRUE,quote = F)
saveRDS(obj, paste0(out_f,"/spatalk.rds"))
############## LR ana ###################
# obj = readRDS('spatalk.rds')
# out_f = './'
r_object = obj@cellpair
df <- data.frame(
  Name = character(),
  cell_sender = character(),
  cell_receiver = character(),
  stringsAsFactors = FALSE
)

for (name in names(r_object)) {
  sender <- r_object[[name]]$cell_sender
  receiver <- r_object[[name]]$cell_receiver
  df <- rbind(df, data.frame(Name = rep(name, length(sender)), cell_sender = sender, cell_receiver = receiver, stringsAsFactors = FALSE))
}
write.csv(df, paste0(out_f,"/cellpair.csv"), row.names = T,quote = F)
write.csv(obj@meta$rawmeta, paste0(out_f,"/spatalk_meta.csv"), row.names = T,quote = F)

# my_dict <- list()
# for (i in seq_along(obj@cellpair)) {
#     # Extract assay name and length
#     assay_name <- names(obj@cellpair)[i]
#     print(assay_name)
#     assay_length <- nrow(obj@cellpair[[i]])
#     print(assay_length)
#     # Extract L1, L2, and L3 labels from assay name
#     labels <- strsplit(assay_name, " -- ")[[1]]
#     sender <- labels[1]
#     rec <- labels[2]
#     if (sender %in% names(my_dict)){
#         my_dict[[sender]] = my_dict[[sender]] + assay_length
#     } else{
#         # initial
#         print('init')
#         my_dict[[sender]] = assay_length
#     }
#     if (rec %in% names(my_dict)){
#         my_dict[[rec]] = my_dict[[rec]] + assay_length
#     } else{
#         # initial
#         print('init')
#         my_dict[[rec]] = assay_length
#     }
#     # break
# }
# my_df <- as.data.frame(my_dict)
# write.table(my_df, file = paste0(out_f,"/CCI_nn_tps.tsv"), sep='\t', quote = F)

# names_list = colnames(my_df)

# df <- matrix(0, nrow = length(names_list), ncol = length(names_list),
#              dimnames = list(names_list, names_list))
# # Convert the matrix to a data frame
# df <- as.data.frame(df)

# for (i in seq_along(obj@cellpair)) {
#     # Extract assay name and length
#     assay_name <- names(obj@cellpair)[i]
#     print(assay_name)
#     assay_length <- nrow(obj@cellpair[[i]])
#     print(assay_length)
#     # Extract L1, L2, and L3 labels from assay name
#     labels <- strsplit(assay_name, " -- ")[[1]]
#     sender <- labels[1]
#     rec <- labels[2]
#     df[sender,rec] = assay_length
#     }
# write.table(df, file = paste0(out_f,"/CCI_n_df.tsv"), sep='\t', quote = F)

# # 压扁 每个rec而言sender的占比 colsum = 1
# colsums <- colSums(df)
# df_send <- t(apply(df, 1, function(x) x / colsums))
# df_send[is.na(df_send)] <- 0
# # e.g. df_send最后一列是Ex5b作为rec，Ex4只有Ex5b一个rec
# # 每个sender而言rec的占比，rowsum = 1
# rowsums <- rowSums(df)
# df_rec <- apply(df, 2, function(x) x / rowsums)
# df_rec[is.na(df_rec)] <- 0
# # e.g. df_rec最后一行是Ex5b作为sender，没有rec细胞，所以全是0
# write.table(df_send, file = paste0(out_f,"/col4Rec_nn_tps.tsv"), sep='\t', quote = F)
# write.table(df_rec, file = paste0(out_f,"/row4Send_nn_tps.tsv"), sep='\t', quote = F)
