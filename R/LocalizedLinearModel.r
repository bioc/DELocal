

# Returns count data keeping the expression columns at the begining of the data frame.
#' Title
#'
#' @param ddsHTSeq
#' @param isNormalized
#' @param exprsnColPattern
#' @param geneAnnot
#'
#' @return
#' @export
#'
#' @examples
getCountsFromDESeqData <-function(ddsHTSeq,isNormalized,exprsnColPattern,geneAnnot){
  ddsHTSeq_count<-DESeq2::estimateSizeFactors(ddsHTSeq)
  ddsHTSeq_count<-as.data.frame( DESeq2::counts(ddsHTSeq_count,normalized=isNormalized))
  ddsHTSeq_count$ensembl_gene_id <- rownames(ddsHTSeq_count)
  ddsHTSeq_count <- merge(ddsHTSeq_count,geneAnnot,by="ensembl_gene_id")
  ddsHTSeq_count[, c(grep(exprsnColPattern, colnames(ddsHTSeq_count)) ,
                     grep(exprsnColPattern, colnames(ddsHTSeq_count), invert = TRUE))]
}

#' Title
#'
#' @param residuals
#' @param selectA
#' @param selectB
#' @param pvalue.fdr
#' @param adjustPara
#' @param lfc
#'
#' @return
#' @export
#'
#' @examples
newMakeContrastBvsA<-function(residuals, selectA, selectB, pvalue.fdr,adjustPara,lfc)
{
  require(limma)
  select<-c(selectA,selectB)
  designA<-c(rep(1,length(selectA)),rep(0,length(selectB)))
  designB<-c(rep(0,length(selectA)),rep(1,length(selectB)))
  design<-cbind(A=designA,B=designB)
  fit<-lmFit(residuals[,select],design)
  cont.mat<-makeContrasts(AvsB=B-A,levels=design)
  fit2<-contrasts.fit(fit,cont.mat)
  fit3<-eBayes(fit2)
  top<-topTable(fit3,adjust=adjustPara,number=nrow(residuals),lfc = lfc)
  pos<-which(top[,"adj.P.Val"]<=pvalue.fdr)
  contrastBvsA<-top[pos,]
  return(contrastBvsA)
}

#' Finds differentially expressed genes by comparing neighboring genes.
#'
#' @param nearest_neighbours
#' @param smrExpt SummarizedExperiment object
#' @param contrast
#'
#' @return
#' @export
#'
#' @examples
DELocal<-
    function(smrExpt,contrast,nearest_neighbours,pDesign,pValue_cut,logFold_cut){

      smrExpt <- DESeq2::DESeqDataSet(smrExpt, design = pDesign)
      smrExpt <- DESeq2::estimateSizeFactors(smrExpt)
      SummarizedExperiment::assays(smrExpt)$normalized_counts = DESeq2::counts(smrExpt, normalized = TRUE)
      exp_mat <- as.data.frame(SummarizedExperiment::assays(smrExpt)$normalized_counts)
      sample_names <- colnames(exp_mat)

      exp_mat$ensembl_gene_id = rownames(exp_mat)

      require(dplyr)
      linear_model <- LocalizedLinearModel(
        exp_mat %>% left_join(SummarizedExperiment::rowData(smrExpt) %>% as.data.frame(), by=c("ensembl_gene_id" ="ensembl_gene_id")),
        sample_names,
        nearest_neighbours+1)

      p_selectA <- which(SummarizedExperiment::colData(smrExpt)[,contrast[1]]== contrast[2])
      p_selectB <- which(SummarizedExperiment::colData(smrExpt)[,contrast[1]]== contrast[3])

      DELocal_table <-
        newMakeContrastBvsA(
          linear_model,
          selectA = p_selectA,
          selectB = p_selectB,
          pvalue.fdr = pValue_cut,
          adjustPara = "BH",
          lfc = logFold_cut
        )

        return(DELocal_table)
    }

LocalizedLinearModel<-
  function(gene_xprsn_annotation,sample_names,nearest_neighbours){
    require(gtools)


    gene_xprsn_annotation <- gene_xprsn_annotation[order(gene_xprsn_annotation$chromosome_name,gene_xprsn_annotation$start_position),]

    linear_models <-  data.frame(matrix(vector(), 0,length(sample_names), dimnames=list(c(), sample_names)), stringsAsFactors=F)

    for (i in 1:nrow(gene_xprsn_annotation)){
      current_chromosom <- gene_xprsn_annotation[i,]$chromosome_name
      current_loc <- gene_xprsn_annotation[i,]$start_position
      current_ensembl_gene_id <- gene_xprsn_annotation[i,]$ensembl_gene_id

      neighbors <- gene_xprsn_annotation[ gene_xprsn_annotation$chromosome_name ==current_chromosom
                                          & abs(gene_xprsn_annotation$start_position - current_loc)< 500000,]

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
          print(paste(i,current_ensembl_gene_id,sep=" "));
        })

      }else {
        linear_models[i,] <- rep(NA,length(sample_names))
      }
      rownames(linear_models)[i] <- as.character(current_ensembl_gene_id)
    }
    return(linear_models)
  }


# plotting linear model of genes local
# plot_LOCAL_lm("ENSMUSG00000095105",anno_trFitENSG_prot[!duplicated(anno_trFitENSG_prot$ensembl_gene_id), ],c(E13_Tooth_Br,E14_Tooth_Br,E14_Jaw_NA),6)
# Setting ylim in R's plot.lm residual plot
# plotlm(heatr1_lm,pch=16,col=c(rep("red",5),rep("green",5),rep("blue",5)),which=1)
plot_LOCAL_lm<-function(pEnsembl_gene_id, gene_xprsn_annotation,sample_names,nearest_neighbours,isNeigbor){
  require(gtools)
  # parameters :: nearest_neighbours,
  gene_xprsn_annotation <- gene_xprsn_annotation[order(gene_xprsn_annotation$chromosome_name,gene_xprsn_annotation$start_position),]

  i <- which(pEnsembl_gene_id==gene_xprsn_annotation$ensembl_gene_id)
  current_chromosom <- gene_xprsn_annotation[i,]$chromosome_name
  current_loc <- gene_xprsn_annotation[i,]$start_position
  current_ensembl_gene_id <- gene_xprsn_annotation[i,]$ensembl_gene_id

  neighbors <- gene_xprsn_annotation[ gene_xprsn_annotation$chromosome_name ==current_chromosom
                                      & abs(gene_xprsn_annotation$start_position - current_loc)< 500000,]

  neighbors$distance <- neighbors$start_position - gene_xprsn_annotation[i,]$start_position
  neighbors <- neighbors[order(abs(neighbors$distance)),]
  neighbors <- na.omit(neighbors)

  if(dim(neighbors)[1] >= 3){
    num_neighbours <- nearest_neighbours
    max_neighbours <- dim(neighbors)[1]
    if(max_neighbours < nearest_neighbours){
      num_neighbours  <- max_neighbours
    }

    y <- as.numeric( log2(neighbors[ 1,c(sample_names)]))
    median_sample_expression <- as.numeric (apply(log2 (neighbors[2:num_neighbours,sample_names]),2,median))
    x <- (log2(neighbors[ 2:num_neighbours,c(sample_names)]))
    sd_local<-apply(x,1,sd)
    mean_local<-apply(x,1,mean)
    CVsquare_local <- (sd_local/mean_local)^2
    var_local<-apply(x,1,var)
    stable_expression <- as.numeric(x[which.min(var_local),])
    if(isNeigbor==TRUE){

      return(data.frame(id=pEnsembl_gene_id,mean_median_gene=mean(median_sample_expression), CV_median_gene= (sd(median_sample_expression)/mean(y)^2) , mean_stable_neighbour=  mean_local[which(CVsquare_local==min(CVsquare_local))] , CV_stable_neighbour=min(CVsquare_local) ))
    }
    tryCatch(  {
      par(mfrow=c(1,1))
      plot(stable_expression,y,col=c(rep("red",5),rep("green",5),rep("blue",3)),main="")#xlim=c(5,8),ylim=c(4.5,8),
      abline(lm( y ~ stable_expression))
      legend(7,5,bty='n', c('E13 Tooth','E14 Tooth','E14 Jaw'),  pt.bg='white', title="",pt.cex=c(1),cex=0.5,col = c("red","green","blue"),pch=1)
      return( lm( y ~ stable_expression))
   }, error= function(err){
      print(paste(i,current_ensembl_gene_id,sep=" "));
    })

  }
  return(data.frame(id=pEnsembl_gene_id,mean_median_gene=FALSE, CV_median_gene= FALSE , mean_stable_neighbour=FALSE , CV_stable_neighbour=FALSE ))
}



#' Title
#'
#' @param exprsn_contr_column
#' @param pExprsn_Location a dataframe where starting columns for expression values and finishing columns correspond to gene_id, chromosome and location
#' @param lfc
#' @param p_value
#' @param gene_reference should be a charactor vector
#'
#' @return
#' @export
#'
#' @examples
optimize_Local <-
  function(exprsn_contr_column,
           pExprsn_Location,
           lfc,
           p_value,
           true_gene_list) {
    results <- list()
    linear_models_list <- list()
    new_tune_neighbour <-
      data.frame(
        neighbours = numeric(),
        Number_of_genes =  numeric(),
        class  = character()
      )
    DE_local_results <-
      data.frame(neighbours = numeric(),
                 method = character(),
                 AUC = numeric())
    neighbor_to_evaluate = 2:15
    l <- for (i in neighbor_to_evaluate)  {
      linear_models_list[[i]] <-
        LocalizedLinearModel(pExprsn_Location[!duplicated(pExprsn_Location$ensembl_gene_id),], colnames(pExprsn_Location)[exprsn_contr_column$exprsn_column], i)
      DELocal_table <-
        newMakeContrastBvsA(
          linear_models_list[[i]],
          selectA = exprsn_contr_column$contst_selectA,
          selectB = exprsn_contr_column$contst_selectB,
          pvalue.fdr = p_value,
          adjustPara = "BH",
          lfc = lfc
        )

      DELocal_table$true_genes <- rownames(DELocal_table) %in% true_gene_list

      pred_DELocal <-
        ROCR::prediction(abs(DELocal_table$logFC), DELocal_table$true_genes)
      DE_local_results <- rbind( DE_local_results,
                                 data.frame(
                                   neighbours = i - 1,
                                   method = "true_genes",
                                   AUC = ROCR::performance(pred_DELocal, "auc")@y.values[[1]]
                                   )
      )
      }


    gene_reference <- dplyr::select(pExprsn_Location, ensembl_gene_id)

    gene_reference$true_genes <- gene_reference$ensembl_gene_id %in% true_gene_list

    # declare this function first
    performance_neighbor <-
      performance_neighbor(
        linear_models_list = linear_models_list,
        neighbor_nums = neighbor_to_evaluate,
        exprsn_contr_column = exprsn_contr_column,
        p_adjusted =  p_value,
        pLfcCut = lfc,
        reference = gene_reference
      )

    results$linear_models_list <- linear_models_list
    #results$new_tune_neighbour <- new_tune_neighbour
    results$DE_local_results <- DE_local_results
    results$performance_neighbor <- performance_neighbor
    results
  }


performance_neighbor<-function(linear_models_list,neighbor_nums,exprsn_contr_column,p_adjusted,pLfcCut,reference){
  nbor_result<-data.frame()
  for (i in neighbor_nums) {
    tmpResult <- newMakeContrastBvsA(linear_models_list[[i]],
                                     selectA= exprsn_contr_column$contst_selectA,
                                     selectB = exprsn_contr_column$contst_selectB,
                                     pvalue.fdr=p_adjusted, adjustPara = "BH", lfc=pLfcCut)
    tmpResult$ensembl_gene_id <- rownames(tmpResult)
    # tmpResult <- annotate_tooth_genes(tmpResult)
    reference[paste(i)]<- reference$ensembl_gene_id %in% tmpResult$ensembl_gene_id
  }

  performance_nbor <- performance_analysis(paste(neighbor_nums), "true_genes",reference)
  performance_nbor_melt<- reshape2::melt(performance_nbor[,c("neighbour","TPR","PPV","MCC")],id.vars ="neighbour" )
  list(prformnc=performance_nbor,melt_prformnc=performance_nbor_melt)
}

#' Title
#'
#' @param top_genes
#' @param Methods_list
#'
#' @return
#' @export
#'
#' @examples
roc_from_topGenes <- function(top_genes,Methods_list) {
  Method_compare_top <- data.frame(top_genes = numeric(),roc_auc=  numeric(), method  = character())

  t_limma <- prediction(abs(Methods_list$limma[1:top_genes,]$logFC), Methods_list$limma[1:top_genes,]$tooth_genes)
  t.roc.perf = performance(t_limma, measure = "tpr", x.measure = "fpr")
  roc_auc <- performance(t_limma,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "limma"))

  RProd_data_top <- Methods_list$RP[1: top_genes,]
  RProd_data_top$logFC <- ifelse( RProd_data_top[,4] < 1, 1/RProd_data_top[,4] , RProd_data_top[,4])
  RProd_data_top <- RProd_data_top[ with(RProd_data_top , order(-logFC)),]

  pred_RProd <- prediction(1-(RProd_data_top$pfp ),RProd_data_top$tooth_genes)
  roc.perf_RProd = performance(pred_RProd, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_RProd,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "RProd"))

  SAM_top <- Methods_list$SAM[1: top_genes,]
  pred_SAM <- prediction(abs(SAM_top$Score),SAM_top$tooth_genes)
  roc.perf_SAM = performance(pred_SAM, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_SAM,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "SAM"))

  FCros_top_genes <- FCros_e13wt_vs_e14wt [ 1:top_genes,]
  pred_FCros <- prediction(abs(FCros_top_genes$roc_f.value),FCros_top_genes$tooth_genes)
  roc.perf_FCros = performance(pred_FCros, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_FCros,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "FCros"))

  DELocal_table_top <- Methods_list$DELocal[1: top_genes,]
  pred_DELocal <- prediction(abs(DELocal_table_top$logFC), DELocal_table_top$tooth_genes)
  roc.perf_delocal = performance(pred_DELocal, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_DELocal,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "DELocal"))

  # DELocal_table_top <- Methods_list$DELocal_expression[1: top_genes,]
  pred_DELocal <- prediction(abs(DELocal_table_top$logFC), DELocal_table_top$tooth_genes)
  # roc.perf_delocal = performance(pred_DELocal, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_DELocal,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "DELocal_expression"))

  Demi_table_top <- Methods_list$Demi[1: top_genes,]
  pred_demi <- prediction(1-(Demi_table_top$P.value), Demi_table_top$tooth_genes)
  roc.perf_demi = performance(pred_demi, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_demi,"auc")@y.values[[1]]
  Method_compare_top <- rbind(Method_compare_top, data.frame(top_genes = top_genes, roc_auc = roc_auc, method = "Demi"))

  Method_compare_top
}

#' Title
#'
#' @param Methods_list
#'
#' @return
#' @export
#'
#' @examples
compare_methods <- function(Methods_list){
  results <- list()
  Method_compare <- data.frame(roc_auc = numeric(), method  = character())
  pred_limma <- prediction(abs(Methods_list$limma$logFC), Methods_list$limma$tooth_genes)
  results$roc.perf = performance(pred_limma, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_limma,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "limma"))


  pred_DELocal <- prediction(abs(Methods_list$DELocal$logFC), Methods_list$DELocal$tooth_genes)
  results$roc.perf_delocal = performance(pred_DELocal, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_DELocal,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "DELocal"))


  #use logFC and check
  pred_RProd <- prediction(1-Methods_list$RP$pfp,Methods_list$RP$tooth_genes)
  results$roc.perf_RProd = performance(pred_RProd, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_RProd,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "RankProd"))

  pred_SAM <- prediction(abs(Methods_list$SAM[ Methods_list$SAM$Q.value < 0.05,]$Score),Methods_list$SAM[ Methods_list$SAM$Q.value < 0.05,]$tooth_genes)
  results$roc.perf_SAM = performance(pred_SAM, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_SAM,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "SAM"))

  pred_FCros <- prediction(abs(Methods_list$FCros[ Methods_list$FCros$p.value < 0.05,]$roc_f.value),Methods_list$FCros[ Methods_list$FCros$p.value < 0.05,]$tooth_genes)
  results$roc.perf_FCros = performance(pred_FCros, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_FCros,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "FCros"))

  pred_Demi <- prediction((1-Methods_list$Demi$P.value),Methods_list$Demi$tooth_genes)
  results$roc.perf_Demi = performance(pred_Demi, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_Demi,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "Demi"))

  results$Method_compare <- Method_compare
  results
}

#' Title
#'
#' @param Methods_list
#'
#' @return
#' @export
#'
#' @examples
compare_methods_best <- function(Methods_list){
  results <- list()
  Method_compare <- data.frame(roc_auc = numeric(), method  = character())

  pred_limma <- prediction(1-(Methods_list$limma$adj.P.Val), Methods_list$limma$tooth_genes)
  results$roc.perf = performance(pred_limma, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_limma,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "limma"))

  pred_DELocal <- prediction(abs(Methods_list$DELocal$logFC), Methods_list$DELocal$tooth_genes)
  results$roc.perf_delocal = performance(pred_DELocal, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_DELocal,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "DELocal"))

  pred_edgeR <- prediction(1-(Methods_list$edgeR$FDR), Methods_list$edgeR$tooth_genes)
  results$roc.perf_edgeR = performance(pred_edgeR, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_edgeR,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "edgeR"))

  pred_DEseq <- prediction(1-(Methods_list$DEseq$padj), Methods_list$DEseq$tooth_genes)
  results$roc.perf_DEseq = performance(pred_DEseq, measure = "tpr", x.measure = "fpr")
  roc_auc<-performance(pred_DEseq,"auc")@y.values[[1]]

  Method_compare <- rbind(Method_compare, data.frame(roc_auc = roc_auc, method = "DEseq"))


  results$Method_compare <- Method_compare
  results
}

#' Title
#'
#' @param Methods_list
#' @param top_gene
#'
#' @return
#' @export
#'
#' @examples
rnaSeq_rank <- function(Methods_list,top_gene){

  Methods_list$limma <- Methods_list$limma[ order(Methods_list$limma$adj.P.Val),]
  Methods_list$DEseq <- Methods_list$DEseq[ order(Methods_list$DEseq$padj),]
  Methods_list$edgeR <- Methods_list$edgeR[ order(Methods_list$edgeR$FDR),]
  Methods_list$DELocal <- Methods_list$DELocal[ order(-abs(Methods_list$DELocal$logFC)),]

  method_rank <- data.frame (
    DELocal = which( Methods_list$DELocal$tooth_genes == TRUE)[1:top_gene] ,
    limma = which( Methods_list$limma$tooth_genes == TRUE)[1:top_gene],
    DEseq = which( Methods_list$DEseq$tooth_genes == TRUE )[1:top_gene],
    edgeR = which( Methods_list$edgeR$tooth_genes == TRUE)[1:top_gene]
  )
}

#' Title
#'
#' @param Methods
#' @param reference
#'
#' @return
#' @export
#'
#' @examples
combo_results <- function(Methods,reference){
  methods_name <-names(Methods)
  for (i in 1:length(methods_name)){
    reference[methods_name[i]]<- reference$ensembl_gene_id %in% Methods[[methods_name[i]]]$ensembl_gene_id
  }
  reference
}

neighbourAnalysis<-
  function(gene_xprsn_annotation,nearest_neighbours){
    require(gtools)
    # parameters :: nearest_neighbours,
    gene_xprsn_annotation <- gene_xprsn_annotation[order(gene_xprsn_annotation$chromosome_name,gene_xprsn_annotation$start_position),]

    results <-  data.frame(ensembl_gene_id=character(),mgi_symbol=character(),isDE_E13_E14 = logical(), isDE_E14_JAW = logical(),
                           num_neighbours = numeric(), DE_E13_E14_neighbours = numeric(), DE_E14_Jaw_neighbours = numeric(),
                           tooth_genes_ECM= character(), stringsAsFactors=F)

    for (i in 1:nrow(gene_xprsn_annotation)){
      current_chromosom <- gene_xprsn_annotation[i,]$chromosome_name
      current_loc <- gene_xprsn_annotation[i,]$start_position
      current_ensembl_gene_id <- gene_xprsn_annotation[i,]$ensembl_gene_id
      current_mgi_symbol <- gene_xprsn_annotation[i,]$mgi_symbol
      current_isDE_E13_E14 <- gene_xprsn_annotation[i,]$DE_E13_E14
      current_isDE_E14_JAW <- gene_xprsn_annotation[i,]$DE_E14Jaw
      current_tooth_genes <- gene_xprsn_annotation[i,]$tooth_genes_ECM

      neighbors <- gene_xprsn_annotation[ gene_xprsn_annotation$chromosome_name ==current_chromosom
                                          & abs(gene_xprsn_annotation$start_position - current_loc)< 500000,]

      neighbors$distance <- neighbors$start_position - gene_xprsn_annotation[i,]$start_position
      neighbors <- neighbors[order(abs(neighbors$distance)),]
      neighbors <- na.omit(neighbors)


      DE_E13_E14_neighbours <- 0
      DE_E14_Jaw_neighbours <- 0
      num_neighbours <- 0

      if(dim(neighbors)[1] > 1){
        num_neighbours <- nearest_neighbours
        max_neighbours <- dim(neighbors)[1]
        if(max_neighbours < nearest_neighbours){
          num_neighbours  <- max_neighbours
        }
        # num_neighbours <- dim(neighbors)[1]
        neighbors <- neighbors[2:num_neighbours,]
        DE_E13_E14_neighbours <- sum(neighbors$DE_E13_E14)
        DE_E14_Jaw_neighbours <- sum(neighbors$DE_E14Jaw)
      }else {
        num_neighbours <- 1
      }
        results <- rbind(results,data.frame(
          ensembl_gene_id=current_ensembl_gene_id,mgi_symbol=current_mgi_symbol,isDE_E13_E14 = current_isDE_E13_E14,
          isDE_E14_JAW = current_isDE_E14_JAW, num_neighbours = num_neighbours-1, DE_E13_E14_neighbours = DE_E13_E14_neighbours,
          DE_E14_Jaw_neighbours = DE_E14_Jaw_neighbours, tooth_genes_ECM = current_tooth_genes
        ))



    }
    return(results)
  }

#' Title
#'
#' @param gene_xprsn_annotation
#' @param nearest_neighbours
#'
#' @return
#' @export
#'
#' @examples
GO_neighbourAnalysis<-
  function(gene_xprsn_annotation,nearest_neighbours){
    require(gtools)

    results <-
      data.frame(
        ensembl_gene_id = character(),
        closest_neighbours = numeric(),
        window_neighbours = numeric(),
        DE_neighbours = numeric(),
        stringsAsFactors = F
      )
    for (i in 1:nrow(gene_xprsn_annotation)){
      if(gene_xprsn_annotation[i,]$isGO){
            current_chromosom <- gene_xprsn_annotation[i,]$chromosome_name
            current_loc <- gene_xprsn_annotation[i,]$start_position
            current_ensembl_gene_id <- gene_xprsn_annotation[i,]$ensembl_gene_id
            neighbors <- gene_xprsn_annotation[ gene_xprsn_annotation$chromosome_name ==current_chromosom
                                                & abs(gene_xprsn_annotation$start_position - current_loc)< 500000,]
            closest_neighbours <-0 # is closest neighbor from same GO/functional group .
            window_neighbours <-0
            index <- which(neighbors$ensembl_gene_id == current_ensembl_gene_id)
            if(index>1 && neighbors[index-1,]$isGO ){
              closest_neighbours <- closest_neighbours+1
            }
            if(index < dim(neighbors)[1] && neighbors[index+1,]$isGO){
              closest_neighbours <- closest_neighbours+1
            }
            neighbors$distance <- neighbors$start_position - gene_xprsn_annotation[i,]$start_position
            neighbors <- neighbors[order(abs(neighbors$distance)),]
            neighbors <- na.omit(neighbors)

            if(dim(neighbors)[1] > 1){
              num_neighbours <- nearest_neighbours
              max_neighbours <- dim(neighbors)[1]
              if(max_neighbours < nearest_neighbours){
                num_neighbours  <- max_neighbours
              }
              neighbors <- neighbors[2:num_neighbours,]
              window_neighbours <- sum(neighbors$isGO)
              DE_neighbours <- sum(neighbors$DE_gene)

            }else {
              num_neighbours <- 1
            }
            results <-  rbind(results,data.frame(
                  ensembl_gene_id = current_ensembl_gene_id,
                  closest_neighbours = closest_neighbours,
                  window_neighbours = window_neighbours,
                  DE_neighbours = DE_neighbours
                )
              )
        }
      }
    return(results)
  }

getNiegbours <- function(xprsn,gene_name){
  xprsn <- xprsn[order(xprsn$chromosome_name,xprsn$start_position),]
  current_chromosom <- xprsn[ xprsn$mgi_symbol == gene_name,]$chromosome_name
  current_loc <- xprsn[xprsn$mgi_symbol == gene_name,]$start_position
  neighbors <- xprsn[ xprsn$chromosome_name ==current_chromosom & abs(xprsn$start_position -current_loc)< 500000,]
  return(neighbors)
}

#' Title
#'
#' @param Method_name
#' @param reference_name
#' @param reference
#'
#' @return
#' @export
#'
#' @examples
performance_analysis<-function(Method_name, reference_name, reference){
  require(PCRedux)
  results <- data.frame()
  for (i in 1:length(Method_name)){
    x<-performeR( reference[,Method_name[i]],  reference[,reference_name])
    x$neighbour <-Method_name[i]
    results <- rbind(results,x)
  }
  results
}
