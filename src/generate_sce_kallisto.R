#!/usr/bin/env Rscript

suppressPackageStartupMessages( {
  library(SingleCellExperiment)
  library(argparse)
  library(Matrix)
})

parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample - from Kallisto.')

parser$add_argument('--sample',
                    type = "character",
                    help = 'Sample identifier')

parser$add_argument('--working_dir', 
                    type = 'character',
                    help = 'Working directory')

parser$add_argument('--output_fn', 
                    type = 'character',
                    help = 'Output SCE filename (path)')

args <- parser$parse_args()


wd <- args$working_dir
id <- args$sample

# generating sce object

counts <- Matrix::readMM(file.path(wd, 'bustools', id,'output.mtx'))
gene_ids <- readLines(file.path(wd, 'bustools', id,'output.genes.txt'))
barcodes <- readLines(file.path(wd, 'bustools', id,'output.barcodes.txt'))

sce <- SingleCellExperiment(list(counts=t(counts)),
                            colData=DataFrame(Barcode=barcodes),
                            rowData=DataFrame(ID=gene_ids,SYMBOL=gene_ids))
rownames(sce) <- gene_ids


saveRDS(object = sce, file = args$output_fn)
