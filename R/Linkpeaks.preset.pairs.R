## create a list of background peaks  --------------------------------
get_bg_peaks <- function(
  meta.bg = pbmc@assays$peaks2@meta.features,
  my_peak = "chr20-9578684-9578893",
  n_sample = 200) {
  
  gene.chrom <- limma::strsplit2(my_peak,"-")[,1]
  trans.peaks <- row.names(meta.bg)[!grepl(pattern = paste0("^", gene.chrom), x = row.names(meta.bg))]

  MatchRegionStats(meta.feature = meta.bg[trans.peaks,], 
                   query.feature = meta.bg[my_peak,], 
                   features.match = c(
                     "GC.percent",
                     "logcount"
                   ),
                   n = n_sample,
                   verbose = FALSE)
}



## Compute Pearson R and Zscores for a list of links  --------------------------------
Linkpeaks.preset.pairs <- function(
  totest = CRISPR.validation, # has peak and gene_peak columns
  l.bg.peaks = l.bg.peaks,
  my_gene = "PPIF",
  n_sample = 200,
  peak.mat = peak.mat,
  gene.mat = gene.mat) {
  
  # keep all peaks tested for this gene
  my.peaks <- totest %>% dplyr::filter(grepl(gene_peak, pattern = paste0("^",my_gene,"_")))
  gene.all.bg.peaks <- l.bg.peaks[my.peaks$peak]

  # make background peak vector for all tested peaks from that gene
  my_bg_peaks.v <- gene.all.bg.peaks %>% unlist()
  
  #### matrices
  bg_p.counts <- peak.mat[my_bg_peaks.v, ] %>% t()
  gene_counts <- gene.mat[my_gene, ]
  my_p.counts <- peak.mat[my.peaks$peak, ] %>% t()

  #### Pearson correlation test
  if (dim(my_p.counts)[1] == 1) { # dim set for corSparse if only one peak tested
    my_p.counts <- my_p.counts %>% t()
  }

  coef.result <- qlcMatrix::corSparse(X = my_p.counts ,Y = gene_counts %>% Matrix::as.matrix()) %>% as.vector()
  bg.coef <- qlcMatrix::corSparse(X = bg_p.counts ,Y = gene_counts %>% Matrix::as.matrix()) %>% as.data.frame()
  
  bg.coef$gene_peak <- rep(my.peaks$gene_peak, each = n_sample)
  names(coef.result) <- my.peaks$gene_peak
  
  # Z score
  l <- lapply(names(coef.result), function(z){
    bg <- bg.coef %>% dplyr::filter(gene_peak == z) %>% pull(V1)
    (coef.result[z] - mean(x = bg)) / sd(x = bg)
  })
  data.frame(gene_peak = names(unlist(l)),
             gene = my_gene,
             peak = my.peaks$peak,
             Zscore = unlist(l),
             pearson = coef.result
  )
}
