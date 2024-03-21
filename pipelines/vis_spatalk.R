library(SpaTalk)
args <- commandArgs(trailingOnly = TRUE)
# ligand = 'SPP1'
# receptor = 'CD44'
# receptor = 'ITGAV'
# sender = 'malignant'
# receiver = 'Macro'
ligand = args[1]
receptor = args[2]
sender = args[3]
receiver = args[4]
print(receiver)
setwd(args[5])
obj = readRDS('spatalk.rds')
# some name has /, which is not allowed in file name
sender_name <- gsub("/", "_", sender)
receiver_name <- gsub("/", "_", receiver)

pdf(paste0("lrpair_",sender_name,'_',receiver_name,'_',ligand,'_',receptor,".pdf"),width=10,height=8)
out = plot_lrpair(object = obj,ligand = ligand,receptor = receptor,
                    celltype_sender = sender,
                    celltype_receiver = receiver,            
                    size = 4,
                    arrow_length = 0.1)
print(out)
dev.off()

### 4.pathway
pdf(paste0("lr_path_",sender_name,'_',receiver_name,'_',ligand,'_',receptor,".pdf"),width=8,height=8)
out = plot_lr_path(object = obj,ligand = ligand,receptor = receptor,
                   celltype_sender = sender,
                   celltype_receiver = receiver)
print(out)
dev.off()
### 5.path2gene
pdf(paste0("lr_path_gene_",sender_name,'_',receiver_name,'_',ligand,'_',receptor,".pdf"),width=8,height=8)
out = plot_path2gene(object = obj,ligand = ligand,receptor = receptor,
                     celltype_sender = sender,
                     celltype_receiver = receiver)
print(out)
dev.off()