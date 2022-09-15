## SIGNAC Sub-functions for linkpeaks  --------------------------------
LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-")) {
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]
  
  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)
  
  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")
  
  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )
  
  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL
  
  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}

DistanceToTSS <- function(
  peaks,
  genes,
  distance = 200000,
  sep = c("-", "-")
) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}


CollapseToLongestTranscript <- function(ranges) {
  range.df <- data.table::as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}


## SIGNAC LinkPeaks function modified  --------------------------------
# includes multicore and logcounts
require(future)
require(future.apply)
require(qlcMatrix)
require(pbapply)
LinkPeaks2 <- function (object, peak.assay, expression.assay, expression.slot = "data", 
                        gene.coords = NULL, distance = 5e+05, min.distance = NULL, 
                        min.cells = 10, method = "pearson", genes.use = NULL, n_sample = 200, 
                        pvalue_cutoff = 0.05, score_cutoff = 0.05, gene.id = FALSE, 
                        MatchRegionStats.features = c("GC.percent","logcount"),
                        verbose = TRUE, cores_attributed = 1) 
{
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    }
    else if (min.distance == 0) {
      min.distance <- NULL
    }
  }
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = object[[peak.assay]]))
  }
  meta.features <- GetAssayData(object = object, assay = peak.assay, 
                                slot = "meta.features")
  features.match <- c("GC.percent", "count")
  if (!("GC.percent" %in% colnames(x = meta.features))) {
    stop("GC content per peak has not been computed.\n", 
         "Run RegionsStats before calling this function.")
  }
  peak.data <- GetAssayData(object = object, assay = peak.assay, 
                            slot = "counts")
  if (!("count" %in% colnames(x = meta.features))) {
    hvf.info <- FindTopFeatures(object = peak.data)
    hvf.info <- hvf.info[rownames(x = meta.features), ]
    meta.features <- cbind(meta.features, hvf.info)
  }
  expression.data <- GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  peakcounts <- meta.features[rownames(x = peak.data), "count"]
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if (!is.null(x = genes.use)) {
    genes.keep <- intersect(x = names(x = genes.keep[genes.keep]), 
                            y = genes.use)
  }
  expression.data <- expression.data[genes.keep, , drop = FALSE]
  if (verbose) {
    message("Testing ", nrow(x = expression.data), " genes and ", 
            sum(peaks.keep), " peaks")
  }
  genes <- rownames(x = expression.data)
  if (gene.id) {
    gene.coords.use <- gene.coords[gene.coords$gene_id %in% 
                                     genes, ]
    gene.coords.use$gene_name <- gene.coords.use$gene_id
  }
  else {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% 
                                     genes, ]
  }
  if (length(x = gene.coords.use) == 0) { 
    stop("Could not find gene coordinates for requested genes")
  }
  if (length(x = gene.coords.use) < nrow(x = expression.data)) {
    message("Found gene coordinates for ", length(x = gene.coords.use), 
            " genes")
  }
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(peaks = peaks, genes = gene.coords.use, 
                                        distance = distance)
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(peaks = peaks, 
                                              genes = gene.coords.use, distance = min.distance)
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n", "Have you set the proper genome and seqlevelsStyle for ", 
         peak.assay, " assay?")
  }
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  peak.data <- t(x = peak.data)
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  
  
  plan(multiprocess, workers = cores_attributed) # for paralellisation in windows
  
  if (nbrOfWorkers() > 1) { # added package
    mylapply <- future_lapply
  }
  else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  res <- mylapply(X = seq_along(along.with = genes.use), FUN = function(i) {
    peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
    gene.expression <- t(x = expression.data[genes.use[[i]], 
                                             , drop = FALSE])
    gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
    if (sum(peak.use) < 2) {
      return(list(gene = NULL, coef = NULL, zscore = NULL))
    }
    else {
      peak.access <- peak.data[, peak.use, drop = FALSE]
      coef.result <- corSparse(X = peak.access, Y = gene.expression)
      rownames(x = coef.result) <- colnames(x = peak.access)
      coef.result <- coef.result[abs(x = coef.result) > 
                                   score_cutoff, , drop = FALSE]
      if (nrow(x = coef.result) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        peaks.test <- rownames(x = coef.result)
        trans.peaks <- all.peaks[!grepl(pattern = paste0("^", 
                                                         gene.chrom), x = all.peaks)]
        meta.use <- meta.features[trans.peaks, ]
        pk.use <- meta.features[peaks.test, ]
        
        
        
        bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)), 
                           FUN = function(x) {
                             MatchRegionStats(meta.feature = meta.use, 
                                              query.feature = pk.use[x, , drop = FALSE], 
                                              features.match = MatchRegionStats.features,
                                              n = n_sample, verbose = FALSE)
                           })
        bg.access <- peak.data[, unlist(x = bg.peaks), 
                               drop = FALSE]
        bg.coef <- corSparse(X = bg.access, Y = gene.expression)
        rownames(bg.coef) <- colnames(bg.access)
        zscores <- vector(mode = "numeric", length = length(x = peaks.test))
        for (j in seq_along(along.with = peaks.test)) {
          coef.use <- bg.coef[(((j - 1) * n_sample) + 
                                 1):(j * n_sample), ]
          z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
          zscores[[j]] <- z
        }
        names(x = coef.result) <- peaks.test
        names(x = zscores) <- peaks.test
        zscore.vec <- c(zscore.vec, zscores)
        gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
        coef.vec <- c(coef.vec, coef.result)
      }
      gc(verbose = FALSE)
      pval.vec <- pnorm(q = -abs(x = zscore.vec))
      links.keep <- pval.vec < pvalue_cutoff
      if (sum(x = links.keep) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        gene.vec <- gene.vec[links.keep]
        coef.vec <- coef.vec[links.keep]
        zscore.vec <- zscore.vec[links.keep]
        return(list(gene = gene.vec, coef = coef.vec, 
                    zscore = zscore.vec))
      }
    }
  })
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                                3))
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(along.with = unique(x = names(x = coef.vec)))
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = coef.vec)], 
                              x = coef.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  z.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = zscore.vec)], 
                           x = zscore.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  Links(object = object[[peak.assay]]) <- links
  return(object)
}

