#' Finds differentially expressed genes by comparing neighboring genes
#'
#' @param nearest_neighbours How many nearest neighbours within 1 Mb window to evaluate?
#' @param pSmrExpt SummarizedExperiment object
#' @param pDesign design formula
#' @param pValue_cut cut off value for adjusted p-value
#' @param pLogFold_cut cut off value for relative log fold change compared to neighbouring genes
#' @return A data.frame with top significant genes with the following columns:
#'
#' relative.logFC: relative logFC compared to neighbouring genes
#'
#' P.Value: raw p-value
#'
#' adj.P.Value: adjusted p-value
#'
#' B:   log-odds that the gene is differentially expressed
#'
#' @export
#' @importFrom dplyr %>% left_join
#' @importFrom stats model.matrix
#' @importFrom DESeq2 DESeqDataSet estimateSizeFactors counts
#' @importFrom SummarizedExperiment assays rowData colData
#' @examples count_matrix <- as.matrix(read.table(file = system.file("extdata",
#'                                                                   "tooth_RNASeq_counts.txt",
#'                                                                   package = "DELocal")))
#' colData <- data.frame(condition=gsub("\\..*",x=colnames(count_matrix),
#'                                      replacement = ""))
#' gene_location <- read.table(file = system.file("extdata", "gene_location.txt",
#'                                     package = "DELocal"))
#' smrExpt <- SummarizedExperiment::SummarizedExperiment(
#'                                         assays=list(counts=count_matrix),
#'                                         rowData = gene_location,
#'                                         colData=colData)
#' contrast= c("condition","ME13","ME14")
#' require(dplyr)
#' x_genes <- SummarizedExperiment::rowData(smrExpt) %>%
#'       as.data.frame() %>%
#'       filter(chromosome_name=="X") %>% rownames()
#' DELocal_result <- DELocal(pSmrExpt = smrExpt[x_genes,],
#'                          nearest_neighbours = 5, pDesign = ~ condition,
#'                          pValue_cut = 0.05, pLogFold_cut = 0)
DELocal<-
  function(pSmrExpt,nearest_neighbours,pDesign,pValue_cut = 0.05, pLogFold_cut = 0){
    stopifnot("`pSmrExpt` must be a SummarizedExperiment object" = inherits(pSmrExpt,"SummarizedExperiment"))
    if("neighbors_start" %in% (SummarizedExperiment::rowData(pSmrExpt) %>% colnames() )){
        message("User provided neighborhood will be used")
    } else {
        message("Default 1Mb neighborhood will be used")
    }

    pSmrExpt <- DESeqDataSet(pSmrExpt, design = pDesign)
    pSmrExpt <- estimateSizeFactors(pSmrExpt)
    SummarizedExperiment::assays(pSmrExpt)$normalized_counts <- counts(pSmrExpt, normalized = TRUE)
    exp_mat <- as.data.frame(SummarizedExperiment::assays(pSmrExpt)$normalized_counts)
    sample_names <- colnames(exp_mat)

    exp_mat$Xgene_id <- rownames(exp_mat)
    xRowData <- rowData(pSmrExpt) %>% as.data.frame()
    xRowData$Xgene_id <- rownames(xRowData)
    exp_mat <- exp_mat %>% left_join(
        xRowData,
        by = c("Xgene_id" = "Xgene_id")
    )

    linear_model <- .LocalizedLinearModel(
      exp_mat, sample_names,nearest_neighbours + 1
    )

    design_matrix <- model.matrix(pDesign,colData(pSmrExpt))
    DELocal_table <- .DELocal_topTable(pLinear_model = linear_model, pDesign_matrix = design_matrix,
                                      pLogFold_cut = pLogFold_cut,pValue_cut = pValue_cut)
    colnames(DELocal_table)[1] <- "relative.logFC"

    return(DELocal_table[,c("relative.logFC","P.Value","adj.P.Val","B")])
  }

# Table of Top Genes from DELocal Linear Model Fit
#' @import limma
.DELocal_topTable <- function(pLinear_model,pDesign_matrix,pLogFold_cut,pValue_cut) {
    fit <- lmFit(pLinear_model,pDesign_matrix)
    fit2 <- contrasts.fit(fit, coefficients = ncol(pDesign_matrix))
    fit3 <- eBayes(fit2)
    top <- limma::topTable(fit3,adjust="BH",number=nrow(pLinear_model),lfc = pLogFold_cut)
    pos <- which(top[,"adj.P.Val"]<=pValue_cut)
    return(top[pos,])
}


#' Returns median expression from different conditions of genes from a neighbourhood of a gene of interest
#'
#' @param pSmrExpt SummarizedExperiment object
#'
#' @param pNearest_neighbours How many nearest neighbours within 1 Mb window to plot
#' @param pDesign design formula
#' @param colorFactor The coloring factor
#' @param pGene_id The gene of interest
#' @return a list which contains both the data from the neighbourhood and a ggplot object
#' @export
#'
#' @importFrom DESeq2 DESeqDataSet estimateSizeFactors counts
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom dplyr %>% left_join
#' @import ggplot2
#' @examples count_matrix <- as.matrix(read.table(file = system.file("extdata",
#'                                                                   "tooth_RNASeq_counts.txt",
#'                                                                   package = "DELocal")))
#' colData <- data.frame(condition=gsub("\\..*",x=colnames(count_matrix),
#'                                      replacement = ""))
#' gene_location <- read.table(file = system.file("extdata", "gene_location.txt",
#'                                     package = "DELocal"))
#' smrExpt <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=count_matrix),
#'                                             rowData = gene_location,
#'                                             colData = colData)
#' contrast= c("condition","ME13","ME14")
#' require(dplyr)
#' x_genes <- SummarizedExperiment::rowData(smrExpt) %>%
#'       as.data.frame() %>%
#'       filter(chromosome_name=="X") %>% rownames()
#' DELocal::plotNeighbourhood(pSmrExpt = smrExpt, pGene_id = "ENSMUSG00000059401")
plotNeighbourhood<- function(pSmrExpt, pNearest_neighbours=5, pDesign = ~ condition,
                             colorFactor = "condition", pGene_id){
    stopifnot("`pSmrExpt` must be a SummarizedExperiment object" = inherits(pSmrExpt,"SummarizedExperiment"))
    if(!(pGene_id %in% rownames(pSmrExpt))){
        message(paste(pGene_id," does not exist"))
        return()
    }

    pSmrExpt <- DESeqDataSet(pSmrExpt , design = pDesign)
    pSmrExpt <- estimateSizeFactors(pSmrExpt)
    SummarizedExperiment::assays(pSmrExpt)$normalized_counts = counts(pSmrExpt, normalized = TRUE)

    selected_gene <- pSmrExpt[pGene_id] %>% rowData()
    if ("neighbors_start" %in% (rowData(pSmrExpt) %>% colnames())) {
        message("User provided neighborhood will be used")
    } else {
        # message("Default 1Mb neighborhood will be used")
        selected_gene$neighbors_start <- selected_gene$start_position-500000
        selected_gene$neighbors_end <- selected_gene$start_position+500000
    }
    neighborsSmExp <- .getNeighbours(pSmrExpt = pSmrExpt,
                                    pSelected_gene = selected_gene,
                                    pNearest_neighbours = pNearest_neighbours )

    exp_mat <- as.data.frame(SummarizedExperiment::assays(neighborsSmExp)$normalized_counts)
    exp_mat$ensembl_gene_id <- rownames(exp_mat)
    sample_names <- colnames(exp_mat)

    states <- colData(pSmrExpt) %>% as.data.frame() %>% dplyr::pull(colorFactor) %>% unique()

    for(i in states){
        xx <-colData(pSmrExpt) %>% as.data.frame() %>% dplyr::filter(!!sym(colorFactor)==i) %>% rownames()
        exp_mat[,i] <- exp_mat[,xx] %>% as.matrix() %>% matrixStats::rowMedians()
    }

    result_data <- exp_mat[,c("ensembl_gene_id",as.character(states))] %>%
        dplyr::left_join(as.data.frame(rowData(neighborsSmExp)),
                         by = c("ensembl_gene_id" = "ensembl_gene_id"))

    results <- list()
    results$data <-  reshape2::melt(result_data[, c(as.character(states),
                                     c("start_position","ensembl_gene_id", "chromosome_name"))],
                                    id = c("start_position","ensembl_gene_id", "chromosome_name"))

    results$plot <- ggplot(results$data,aes(x = start_position, y = value, text = ensembl_gene_id)) +
        geom_point(aes(colour = variable),size = 3) +
        scale_x_continuous("Gene start distance ",
                           breaks = results$data$start_position,
                           labels = results$data$ensembl_gene_id,
                           limits = c(selected_gene$neighbors_start, selected_gene$neighbors_end)) +
        ylab("Normalized rna-seq") + ggtitle(paste(selected_gene$ensembl_gene_id)) +
        theme_bw() + theme(axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
            legend.position = "bottom",
            legend.title = element_blank())
    message(results$plot)
    return(results)
}

#' @importFrom stats na.omit
#' @importFrom stats lm median
.LocalizedLinearModel<-
  function(gene_xprsn_annotation,sample_names,nearest_neighbours){

    gene_xprsn_annotation <- gene_xprsn_annotation[order(gene_xprsn_annotation$chromosome_name,
                                                         gene_xprsn_annotation$start_position),]

    linear_models <-  data.frame(matrix(vector(), 0,length(sample_names),
                                        dimnames=list(c(), sample_names)), stringsAsFactors=FALSE)

    for (i in seq_len(nrow(gene_xprsn_annotation))){
      current_chromosom <- gene_xprsn_annotation[i,]$chromosome_name
      current_ensembl_gene_id <- gene_xprsn_annotation[i,]$Xgene_id

      if ("neighbors_start" %in% colnames(gene_xprsn_annotation)) {
        current_neighbors_start <- gene_xprsn_annotation[i,]$neighbors_start
        current_neighbors_end <- gene_xprsn_annotation[i,]$neighbors_end
        neighbors <- gene_xprsn_annotation %>% dplyr::filter(
            chromosome_name == current_chromosom &
            start_position >= current_neighbors_start  &
            start_position <= current_neighbors_end
        )
      } else{
        current_loc <- gene_xprsn_annotation[i,]$start_position
        neighbors <- gene_xprsn_annotation[gene_xprsn_annotation$chromosome_name == current_chromosom
                                & abs(gene_xprsn_annotation$start_position - current_loc) < 500000, ]
      }
      neighbors$distance <- neighbors$start_position - gene_xprsn_annotation[i,]$start_position
      neighbors <- neighbors[order(abs(neighbors$distance)),]
      neighbors <- na.omit(neighbors)

      if(dim(neighbors)[1] >= 2){
        num_neighbours <- nearest_neighbours
        max_neighbours <- dim(neighbors)[1]
        if(max_neighbours < nearest_neighbours){
          num_neighbours  <- max_neighbours
        }

        y <- as.numeric( (neighbors[ 1,c(sample_names)]))
        median_sample_expression <- as.numeric (apply( (neighbors[2:num_neighbours,sample_names]),2,median))
        tryCatch(  {lm<-lm( y ~ median_sample_expression)
              linear_models[i,] <- lm$residuals
        }, error= function(err){
          message(paste(i,current_ensembl_gene_id,sep=" "));
        })

      }else {
        linear_models[i,] <- rep(NA,length(sample_names))
      }
      rownames(linear_models)[i] <- as.character(current_ensembl_gene_id)
    }
    return(linear_models)
  }

#' @importFrom stats na.omit
.getNeighbours <- function(pSmrExpt,pSelected_gene,pNearest_neighbours){
  neighborsSmExp <- subset(pSmrExpt, subset = (chromosome_name == pSelected_gene$chromosome_name &
                               start_position >= pSelected_gene$neighbors_start  &
                               start_position <= pSelected_gene$neighbors_end))
  neighbors <- SummarizedExperiment::rowData(neighborsSmExp) %>% as.data.frame()
  neighbors$distance <- neighbors$start_position - pSelected_gene$start_position
  neighbors <- neighbors[order(abs(neighbors$distance)),]
  neighbors <- na.omit(neighbors)
  if(dim(neighborsSmExp)[1] > pNearest_neighbours ){
      neighbors<-neighbors[seq_len(pNearest_neighbours),]
  }
  return(neighborsSmExp[neighbors$ensembl_gene_id,])
}

utils::globalVariables(c("chromosome_name", "start_position", "value",
                         "ensembl_gene_id","variable"))
