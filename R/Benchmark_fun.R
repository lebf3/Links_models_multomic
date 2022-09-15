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


## SIGNAC LinkPeaks function modified : Pearson without null -------------------
# includes multicore and logcounts

object = pbmc
peak.assay = "peaks2"
expression.assay = "RNA" 
expression.slot = "data" 
gene.coords = NULL
distance = 5e+05
min.distance = NULL
min.cells = 10
method = "pearson"
genes.use = NULL
n_sample = 200
score_cutoff = 0.05
gene.id = FALSE
verbose = TRUE
cores_attributed = 1

require(future)
require(future.apply)
require(qlcMatrix)
LinkPeaks.P <- function (object, peak.assay, expression.assay, expression.slot = "data", 
                        gene.coords = NULL, distance = 5e+05, min.distance = NULL, 
                        min.cells = 10, method = "pearson", genes.use = NULL, n_sample = 200, 
                        score_cutoff = 0.05, gene.id = FALSE, 
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
  
  peak.data <- GetAssayData(object = object, assay = peak.assay, 
                            slot = "counts")
  
  expression.data <- GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  
  peakcounts <- rowSums(x = peak.data > 0)
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
    # gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
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
        return(data.frame(gene = genes.use[[i]], 
                          peak = NA, 
                          PearsonR = NA))
      }
    }
    return(data.frame(gene = genes.use[[i]], 
                      peak = row.names(coef.result), 
                      PearsonR = coef.result))
  })
}

