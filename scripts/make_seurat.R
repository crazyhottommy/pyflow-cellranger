library(Seurat)
library(dplyr)
## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


# https://github.com/satijalab/seurat/issues/671
# run NormalizeData first and then FindVariableGenes
PreprocessData<- function(object,
                                x.low.cutoff = 0.05,
                                x.high.cutoff = 10,
                                y.cutoff = 0.5,
                                num.pc = 20,
                                pc.use = NULL,
                                do.par =TRUE,
                                num.cores = 2,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = 0.8,
                                k.param = 30,
                                ...){
		object<- NormalizeData(object = object)

        object<- FindVariableGenes(object = object,
                                   x.low.cutoff = x.low.cutoff,
                                   x.high.cutoff = x.high.cutoff,
                                   y.cutoff = y.cutoff)
        
        object<- ScaleData(object = object, genes.use = object@var.genes,
                           vars.to.regress = c("nUMI"), block.size = 1000,
                           min.cells.to.block=3000,
                           display.progress = TRUE, do.par = TRUE, num.cores = num.cores)

        object<- RunPCA(object = object, pc.genes = object@var.genes,
                        pcs.compute = num.pc, do.print = FALSE)

        if (is.null(pc.use)){
                object<- JackStraw( object = object, num.replicate = 100, num.cores = num.cores,
                                    do.par = T, num.pc = num.pc)

                object <- JackStrawPlot(object = object, PCs = 1:num.pc, score.thresh = score.thresh)

                PC_pvalues<- object@dr$pca@jackstraw@overall.p.values

                ## determin how many PCs to use.
                pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1

        }

        # add significant pc number to metadata, need to have names same as the cells
        pc.use.meta<- rep(pc.use, length(object@cell.names))
        names(pc.use.meta)<- object@cell.names
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object <- RunTSNE(object = object, dims.use = 1:pc.use, do.fast = TRUE)
        object <- FindClusters(object = object, reduction.type = "pca",
                                dims.use = 1:pc.use,
                                n.start = n.start,
                                k.param = k.param,
                                nn.eps = nn.eps, resolution = resolution,
                                print.output = FALSE,
                                save.SNN = TRUE, force.recalc = TRUE)
        return(object)
}


sample<- snakemake@wildcards[["sample"]]
data.dir = snakemake@params[["cellranger_count_dir"]]
# Load the seurat_obj dataset
data <- Read10X(data.dir = data.dir)

seurat_obj <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, 
    project = sample)

make_seurat_obj_args<- snakemake@params[["make_seurat_obj_args"]]
command<- paste("PreprocessData", "(", "seurat_obj,", make_seurat_obj_args, ")")

seurat_obj<- eval(parse(text=command))

saveRDS(seurat_obj, file = paste0("seurat_objs/", sample, "_seurat.rds"))
