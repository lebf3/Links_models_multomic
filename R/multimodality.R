## Returns scaled Pearson R of a link's null distribution  ---------------------
links_distributions <- function(
  my_peak = "chr16-50684843-50685984",
  my_gene = "NOD2",
  n_sample = 1000,
  remove_cor.peaks = F, # removes peaks called in the same cell archetype if specific to it
  RNA.matrix = pbmc[["RNA"]]@data,
  peaks.matrix = pbmc[["peaks2"]]@counts,
  meta = pbmc[["ATACassay"]]@meta.features
) {
  
  if (!remove_cor.peaks) {
    meta.bg <- meta
  }
  else{
    if (meta[my_peak,"peak_call_count"] == 4) { 
      # keep all peaks for non-specific peaks
      meta.bg <- meta
    }
    else{
      # remove columns called in the same cell archetype 
      rem.cols <- meta[my_peak,c("Lymp","NK","Mono","B_cell")] %>% unlist(., use.names=FALSE) 
      rem.cols <- c("Lymp","NK","Mono","B_cell")[rem.cols]
      rem.peaks <- apply(as.data.frame(meta[,rem.cols]) ,1, any)
      meta.bg <- meta[!rem.peaks,]
    }
  }
  
  gene.chrom <- limma::strsplit2(my_peak,"-")[,1]
  trans.peaks <- row.names(meta.bg)[!grepl(pattern = paste0("^", gene.chrom, "-"), x = row.names(meta.bg))]
  
  
  my_bg_peaks <- MatchRegionStats(meta.feature = meta.bg[trans.peaks,], 
                                  query.feature = meta[my_peak,], 
                                  features.match = c(
                                    "GC.percent",
                                    "logcount"
                                  ),
                                  n = n_sample,
                                  verbose = FALSE)
  
  #### matrices
  bg_p.counts <- peaks.matrix[my_bg_peaks, ] %>% t()
  gene_counts <- RNA.matrix[my_gene, ]
  my_p.counts <- peaks.matrix[my_peak,]
  
  #### Pearson correlation test of cis link
  coef.result <- qlcMatrix::corSparse(X = Matrix::as.matrix(my_p.counts),
                                      Y = Matrix::as.matrix(gene_counts)) %>% 
    as.vector()
  
  #### Pearson correlation test of trans links
  bg.coef <- qlcMatrix::corSparse(X = bg_p.counts,Y = Matrix::as.matrix(gene_counts))
  bg.z <- scale(bg.coef)
  
  z <- (coef.result - mean(x = bg.coef)) / sd(x = bg.coef)
  data.frame(bg.z = bg.z, 
             z = z, 
             gene_peak = paste0(my_gene,"_", my_peak),
             bg_peak = my_bg_peaks)
}




## Bimodality test of each link null distributions  ----------------------------
require(multimode)
bimodality_test <- function(
  my_peak = "chr16-50684843-50685984",
  my_gene = "NOD2",
  n_sample = 1000,
  RNA.matrix = pbmc[["RNA"]]@data,
  peaks.matrix = pbmc[["peaks2"]]@counts,
  meta = pbmc[["ATACassay"]]@meta.features
) {
  
  gene.chrom <- limma::strsplit2(my_peak,"-")[,1]
  trans.peaks <- row.names(meta)[!grepl(pattern = paste0("^", gene.chrom, "-"), x = row.names(meta))]
  
  my_bg_peaks <- MatchRegionStats(meta.feature = meta[trans.peaks,], 
                                  query.feature = meta[my_peak,], 
                                  features.match = c(
                                    "GC.percent",
                                    "logcount"
                                  ),
                                  n = n_sample,
                                  verbose = FALSE)
  
  #### matrices
  bg_p.counts <- peaks.matrix[my_bg_peaks, ] %>% t()
  gene_counts <- RNA.matrix[my_gene, ]
  my_p.counts <- peaks.matrix[my_peak,]
  
  #### Pearson correlation test of cis link
  coef.result <- qlcMatrix::corSparse(X = Matrix::as.matrix(my_p.counts),
                                      Y = Matrix::as.matrix(gene_counts)) %>% 
    as.vector()
  
  #### Pearson correlation test of trans links
  bg.coef <- qlcMatrix::corSparse(X = bg_p.counts, 
                                  Y = Matrix::as.matrix(gene_counts))
  bg.z <- scale(bg.coef)
  
  z <- (coef.result - mean(x = bg.coef)) / sd(x = bg.coef)
  
  
  #### testing for multiple modes
  m.t <- modetest(bg.z)
  data.frame(gene_peak = paste0(my_gene, "_", my_peak),
             H0.p = m.t$p.value)
}

