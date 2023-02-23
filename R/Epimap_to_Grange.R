## Liftover Epimap cCRE predictions found in all replicates  --------------------------------
# Download and select Epimap links present in all replicates for this celltype
# wget("https://personal.broadinstitute.org/cboix/epimap/metadata/Imputation_Metadata.xlsx")

Epimap_to_Grange <- function(epimap_celltype = "CD14 MONOCYTE", epi.meta = Epi.meta) {
  require(AnnotationHub)
  
  # query samples for this celltype from meta data
  epi.samples <- epi.meta$BSSID[which(epi.meta$`Extended Info` == epimap_celltype)]
  
  l <- lapply(epi.samples, function(x){
    u <- url(paste0("https://personal.broadinstitute.org/cboix/epimap/links/links_corr_only/",x,"_collated_pred.tsv.gz"))
    txt <- readLines(gzcon(u))
    epi.sample <- read.table(textConnection(txt)) 
    colnames(epi.sample) <- c("chr","start","end","gene","score","state")
    epi.sample <- epi.sample %>% unite(col = peak, chr, start, end, sep = "-", remove = F)
    epi.sample$sample <- x
    epi.sample <- epi.sample %>% unite(col = gene_peak, gene, peak, sep = "_", remove = F)
    epi.sample <- epi.sample[!duplicated(epi.sample$gene_peak),]
    epi.sample
  })
  
  epi.links <- do.call(rbind,l)
  
  # keep gene peak links present in all replicates
  keep.gene_peak <- epi.links %>% 
    group_by(gene_peak) %>% 
    dplyr::count() %>% 
    dplyr::filter(n == length(epi.samples)) %>% 
    pull(gene_peak)
  
  epi.links.filtered <- epi.links %>%
    dplyr::filter(gene_peak %in% keep.gene_peak) %>% 
    distinct(gene_peak,.keep_all = T) %>% 
    dplyr::select(gene_peak:gene)
  
  # get gene names
  ensembl.hsa <- biomaRt::useEnsembl(biomart="ensembl",
                                    dataset="hsapiens_gene_ensembl")
  
  gene.names <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                              "external_gene_name"),
                               filters = "ensembl_gene_id",
                               values = epi.links.filtered$gene ,
                               ensembl.hsa)
  
  epi.links.filtered$gene_name <- gene.names$external_gene_name[match(
    epi.links.filtered$gene, 
    gene.names$ensembl_gene_id
  )]
  
  #### lisftover to hg38
  ahub <- AnnotationHub()
  ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
  chain <- ahub.chain[ahub.chain$title == "hg19ToHg38.over.chain.gz"]
  chain <- chain[[1]]
  
  grl2.hg38 <- liftOver(makeGRangesFromDataFrame(epi.links.filtered, keep.extra.columns = T), chain)
  grl2.hg38 <- unlist(grl2.hg38)
  grl2.hg38$peak_hg38 <- GRangesToString(grl2.hg38)
  grl2.hg38
}