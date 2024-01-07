### ---------------
###
### Create: qyl
### Date: 2023-09-04 16:40:00
### Email: 524583919@qq.com
### Function：
### - the pipeline of WGCNA
### -
###
### Update Log: 2023-09-04
###
### ---------------

#### ----- Function ------ ####
#' the pipeline of WGCNA
#'
#' @title the pipeline of WGCNA
#' @description the pipeline of WGCNA
#' @param exp the input expression matrix
#' @param clinic the input clinic
#' @param num_disc the number of discrete variables
#' @param strict whether to strict
#' @param threads the number of threads
#' @param rate_sap the rate of sample
#' @param num_genes the number of genes
#' @param minModule the number of minModule
#' @param CutHeight the height of cut
#' @param deepSplit the deep of split
#' @param grp_nm the group name
#' @param dir_nm the dir name
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{WGCNA_pipeline}}, \code{\link{miRNA_mRNA}}, \code{\link{Mfuzz_pipeline}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/WGCNA_pipeline.R")
#' WGCNA_pipeline(exp = NULL, # 输入原始的TPM即可，行为基因，列为样本
#'               clinic = NULL,  # 行为样本，列为变量，离散变量在前，连续变量在后
#'               num_disc = NA, strict = T,threads = 8, # 一般默认值就可以
#'               rate_sap = 0.3, num_genes = 5000, # 一般默认值就可以
#'               minModule=30,CutHeight=0.25,deepSplit = 2,
#'               grp_nm = "WGCNA_try1",dir_nm = "WGCNA")
#'               }
#'
#' @import WGCNA
#' @importFrom data.table fread fwrite
#' @importFrom stringr str_replace_all
#' @importFrom stats dist hclust mad
#' @importFrom graphics abline barplot hist par text
#' @export
#'
WGCNA_pipeline <- function(exp = NULL, # 输入原始的TPM即可，行为基因，列为样本
                           clinic = NULL,  # 行为样本，列为变量，离散变量在前，连续变量在后
                           num_disc = NA, strict = T,threads = 8, # 一般默认值就可以
                           rate_sap = 0.3, num_genes = 5000, # 一般默认值就可以
                           minModule=30,CutHeight=0.25,deepSplit = 2,
                           grp_nm = "WGCNA_try1",dir_nm = "WGCNA"){

  ### 1.library ####
  # suppressMessages({
  #   library(magrittr)
  #   library(WGCNA)
  #   library(tidyverse)
  #   library(data.table)
  #   library(stringr)
  # })
  if(strict == T){
    options(stringsAsFactors = FALSE)
    allowWGCNAThreads(nThreads = threads) # 允许R最大线程运行
    enableWGCNAThreads(nThreads = threads) # 打开多线程
  }
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm)
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  unlink(output_dir,recursive = T)
  unlink(photo_dir,recursive = T)
  dir.create(output_dir,recursive = T)
  dir.create(photo_dir,recursive = T)

  ### 2.df_data ####
  df_data <- exp %>%
    na.omit()
  if(strict==T){
    df_data <- df_data %>%
      add(1) %>%
      log2()
  }
  if(strict == T){
    df_data <- df_data[rowSums(df_data)>rate_sap*ncol(df_data),]
    df_data <- df_data[order(apply(df_data,1,mad), decreasing = T)[1:num_genes],]
  }

  ### 3.检查缺失、离群样本 ####
  ## 3.1 缺失样本
  if(strict == T){
    gsg <-  goodSamplesGenes(df_data, verbose = 3) # 检测缺失值
    gsg$allOK
    if (!gsg$allOK) { # 为F，进行QC
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
        print(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
      if (sum(!gsg$goodSamples)>0)
        print(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
      # Remove the offending genes and samples from the data:
      datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
    }
  }

  ## 3.2 离群样本
  par(mfrow=c(1,1), # 仅仅只画出聚类图
      pin=c(10,8),
      mai=c(1,0.5,1,0.5))
  list_hclust <- df_data %>% # Dendrogram
    t() %>%
    dist(method = "euclidean") %>%
    hclust(method = "average")
  pdf(file = paste0(photo_dir,'/1_Sample_clustering_to_detect_outliers.pdf'),
      width = 10, height = 6)
  plot(list_hclust,
       main = "Sample clustering to detect outliers", sub="", xlab="",
       hang = -1,
       cex.lab = 1, cex.axis = 1, cex.main = 1.5,cex=0.4)
  abline(h = 20000, col = "red")
  dev.off()

  df_fil <- df_data # 也可以不过滤
  datExpr <- t(df_fil) # !!! 行为样本

  ## 3.3 clinic_data
  df_group <- clinic
  ord_group <- match(colnames(df_fil),
                     rownames(df_group)) # 将df_group,按照df_fil排序
  df_group <- df_group[ord_group,,drop = FALSE]
  # df_group[is.na(df_group)] <- 0

  ### 4.选取软阈值K ####
  ## 4.1 最佳软阈值计算
  powers <-  c(1:10, seq(from = 12, to=20, by=2))
  sft <-  pickSoftThreshold(datExpr, # 行为样本
                            powerVector = powers,
                            verbose = 5)
  sft$fitIndices
  sft$powerEstimate

  pdf(file = paste0(photo_dir,'/2_Soft_Threshold.pdf'),
      width = 7,height = 5)
  {
    par(mfrow = c(1,2))
    plot(x = sft$fitIndices[,1],# power_value
         y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = "Soft Threshold (power_value)",
         ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels = powers,
         cex = 0.9,
         col = "red")
    abline(h = 0.8, col = "red") # 一般使用0.8
    plot(sft$fitIndices[,1], sft$fitIndices[,5], # mean.k，相同R^2的情况下，mean.k越大越好
         xlab="Soft Threshold (power_value)",
         ylab="Mean Connectivity",
         type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5],
         labels=powers,
         cex=0.9,
         col="red")
  }
  dev.off()

  ## 4.2 软阈值的手动选择
  power_value <- sft$powerEstimate # 使用推荐的值
  if (is.na(power_value)){
    cat("power_mode: NA\n",sep = "",
        file = paste0(output_dir,"/WGCNA_basic_info.txt"))

    nSamples <- ncol(df_data)
    type <- "unsigned"
    power_value <-  ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                           ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                                  ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                         ifelse(type == "unsigned", 6, 12))
                           )
    )
  }else{
    cat("power_mode: powerEstimate\n",sep = "",
        file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  }
  cat("power_value: ",power_value,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  k <- softConnectivity(datExpr=datExpr,
                        power=power_value)
  pdf(file = paste0(photo_dir,"/3_scale_free.pdf"),
      width =8,height = 5 )
  {
    # 指定β下网络是否scale free
    par(mfrow=c(1,2))
    hist(k) # 无尺度网络的条形图分布
    scaleFreePlot(k,main="Check Scale free topology\n")
  }
  dev.off()

  ### 5.网络模块构建 ####
  ## 5.1 模块构建
  cor <-  WGCNA::cor # 必须，blockwiseModules会报错
  net <-  blockwiseModules(
    datExpr, # 列为基因
    power=power_value,
    maxBlockSize = num_genes, # 最大的block的大小，最好比所有基因大，保证分到一个block
    TOMType = "unsigned",
    minModuleSize = minModule, # 最小模块的基因数，值越小保留小模块，default=30
    mergeCutHeight = CutHeight, # 合并相似性模块的距离，值越小不容易合并，模块多，defaul=0.25
    deepSplit = deepSplit, # 模块检测时进行递归分裂的深度,划分模块的敏感度，值越大越敏感，模块越多,defaul=2
    numericLabels = TRUE, # 后面可以再转换为颜色
    pamRespectsDendro = FALSE,
    saveTOMs = F, # 作用不大，不保存了，大概增加15s
    saveTOMFileBase = paste0(output_dir,"/TOM"), # 存储文件名

    ##
    nThreads=threads,
    verbose = 3
  )
  cor <-  stats::cor
  save(net,file = paste0(output_dir,"/net.RData"))
  cat("strict: ",strict,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("rate_sap: ",rate_sap,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("minModuleSize: ",minModule,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("mergeCutHeight: ",CutHeight,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("deepSplit: ",deepSplit,"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))

  ## 5.2 模块信息
  # 转换标签为颜色绘制
  color_vec <- labels2colors(sort(unique(net$colors)))  # 将颜色数字0-n转为标签展示
  moduleColors <-  labels2colors(net$colors)
  table(moduleColors)[color_vec] # 查看所有模块里基因的数量
  cat("gene_number: ",length(net$colors),"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("merged_module_num: ",length(table(net$colors)),"\n",sep = "",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("unmerged_module:",table(net$unmergedColors),"\n",sep = " ",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("merged_module:",table(net$colors),"\n",sep = " ",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))
  cat("merged_module_color:",color_vec,"\n",sep = " ",append = T,
      file = paste0(output_dir,"/WGCNA_basic_info.txt"))

  module_num <- data.frame(module=color_vec,num=table(net$colors))
  write.table(module_num,paste0(output_dir,"/module_num.txt"),
              sep = "\t",row.names = F,col.names = T,quote = F)

  ## 5.3 模块树状图
  for (i in 1:length(net$dendrograms)) {
    pdf(file = paste0(photo_dir,'/4_cluster_dendrogram.pdf'),
        width=8,height=5)
    print(
      plotDendroAndColors(
        net$dendrograms[[i]],
        moduleColors[net$blockGenes[[i]]],
        "Module colors", # 绘制树状图和下面的模块颜色
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05
      )
    )
    dev.off()
  }

  ### 6.模块间相关性 ####
  ## 6.1 模块相似性矩阵
  MEs_col <- MEs <-  net$MEs
  colnames(MEs_col) <-  paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col <-  orderMEs(MEs_col)
  colnames(MEs)
  colnames(MEs_col)

  ## 6.2 模块相似性树状图和热图
  pdf(file = paste0(photo_dir,"/5_eigengene_adjacency_heatmap.pdf"),
      width = 6,height = 6)
  plotEigengeneNetworks(
    MEs_col,
    "Eigengene adjacency heatmap",
    marDendro = c(3,3,2,4),
    marHeatmap = c(3,4,2,2),
    plotDendrograms = T,
    xLabelsAngle = 90
  )
  dev.off()

  ### 7.模块与性状关系 ####
  ## 7.1 计算MEs与clinical相关性
  {
    # 离散变量
    design <- as.data.frame(matrix(numeric(0),nrow = nrow(df_group)))
    if(num_disc!= 0 & !is.na(num_disc) & !is.null(num_disc)){
      for(i in 1:num_disc){
        grp_disc <- df_group[,i]
        i_design <- model.matrix(~0+ grp_disc)
        colnames(i_design) <-  grp_disc %>%
          as.factor() %>%
          levels()
        design <- cbind(design,i_design) # 构建一列或者多列的设计矩阵
      }
    }

    # 连续变量
    # 0,1表示的二分类变量也可以当作连续变量，或者哑变量化的
    if(ncol(df_group)>1){
      design <- cbind(design,df_group[(num_disc+1):ncol(df_group)])
    }
  }
  {
    # MEs矩阵
    MEs <- MEs_col
  }
  {
    # 计算排序后MEs与clinical的相关性
    moduleTraitCor <-  stats::cor(
      MEs,
      design , # 每列可以是数值型的，或者是哑变量化的
      use = "pairwise.complete.obs",
      method = "spearman"
    )
    moduleTraitPvalue <-  corPvalueStudent(
      moduleTraitCor,  # 对moduleTraitCor添加P值
      nrow(datExpr)
    )
    {
      # 生成textMatrix，绘图用显示值用
      textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "")
      dim(textMatrix) <-  dim(moduleTraitCor)
      colnames(textMatrix) <- colnames(moduleTraitCor)
      rownames(textMatrix) <- rownames(moduleTraitCor)
      textMatrix <- as.data.frame(textMatrix)
    }
    {
      # 生成textMatrix1，保存用
      textMatrix1 <-  paste(signif(moduleTraitCor, 2), "(",
                            signif(moduleTraitPvalue, 1), ")", sep = "")
      dim(textMatrix1) <-  dim(moduleTraitCor)
      colnames(textMatrix1) <- colnames(moduleTraitCor)
      rownames(textMatrix1) <- rownames(moduleTraitCor)
      textMatrix1 <- as.data.frame(textMatrix1)
    }
  }
  write.table(moduleTraitCor,file = paste0(output_dir,"/moduleTraitCor.txt"),
              row.names = T,col.names = NA,sep="\t",quote = F)
  write.table(moduleTraitPvalue,file = paste0(output_dir,"/moduleTraitPvalue.txt"),
              row.names = T,col.names = NA,sep="\t",quote = F)
  write.table(textMatrix1,file = paste0(output_dir,"/moduleTrait_textMatrix1.txt"),
              row.names = T,col.names = NA,sep="\t",quote = F)

  ## 7.2 绘制模块与性状相关性图
  pdf(file = paste0(photo_dir,"/6_module_trait_relationships_all.pdf"),
      width = 3+ncol(moduleTraitCor)*0.65,height = 2+nrow(moduleTraitCor)*1.8)
  par(mfrow=c(1,1), mar=c(7, 10, 2, 1))
  labeledHeatmap(Matrix = moduleTraitCor, # 在热图中显示模块和条件之间的相关性
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50), # greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()

  ## 7.3 all module的平均基因表达
  moduleTraitCor
  dir.create(paste0(output_dir,"/module_gene"))
  dir.create(paste0(photo_dir,"/module_gene"))
  for(i_color in unique(moduleColors)){
    module <- i_color # 选择与所研究条件相关性高的module
    dat <- datExpr[,moduleColors==module] %>%
      scale() %>%
      t()

    pdf(file = paste0(photo_dir,"/module_gene/",module,"_eigengene_expression.pdf"),
        width = 12,height = 6)
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
    plotMat(dat,
            nrgcols=30,
            rlabels=F,
            rcols=module,
            main=module,
            cex.main=2)
    barplot(MEs[,paste0("ME",module)],
            col=module,cex.main=2,border = NA,
            main="",ylab="eigengene expression",xlab="samples")
    dev.off()

    module_gene <- colnames(datExpr)[moduleColors==module] # 每一个模块里的基因
    write.table(module_gene,
                file = paste0(output_dir,"/module_gene/",module,"_gene.csv"),
                quote =F,row.names = F,col.names = F)
  }

}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/WGCNA_pipeline.R")
  WGCNA_pipeline(exp = NULL, # 输入原始的TPM即可，行为基因，列为样本
                 clinic = NULL,  # 行为样本，列为变量，离散变量在前，连续变量在后
                 num_disc = NA, strict = T,threads = 8,
                 rate_sap = 0.3, num_genes = 5000, # 一般默认值就可以
                 minModule=30,CutHeight=0.25,deepSplit = 2,
                 grp_nm = "WGCNA_try1",dir_nm = "WGCNA")
}
