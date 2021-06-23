# run iteratively for each dataset that will be integrated
# uses the "process_gene_names.R" function to correct gene names
# and saves each dataset in a new file
# takes as arguments each dataset path and its abbreviated name

source("utils/process_gene_names.R")

args = commandArgs(TRUE)
data = args[1]
name = args[2]

print(name)

data <- readRDS(data)

# make sure we only have the RNA assay
DefaultAssay(data) <- "RNA"
if ( length(Assays(data)) > 1) {data[["SCT"]] <- NULL}

print(data)

processed <- process_gene_names(data)

print(processed)

saveRDS(processed, paste0("output/", name, "_filtered_genes.rds"))
