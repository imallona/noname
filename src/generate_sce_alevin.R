#!/usr/bin/env R

suppressPackageStartupMessages( {
  library(SingleCellExperiment)
  library(argparse)
  library(Matrix)
  library(tximeta) # need to add to environment for R 
})


parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample - from alevin.')

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

# read count matrix
dir<-file.path(wd, 'alevin', id)
files<-file.path(dir, "alevin", "quants_mat.gz")
file.exists(files)

se<- tximeta(files, type="alevin", alevinArgs=list(filterBarcodes=TRUE),txOut = TRUE,skipMeta = TRUE)

sce <- as(se, "SingleCellExperiment")

saveRDS(object = sce, file = args$output_fn)
