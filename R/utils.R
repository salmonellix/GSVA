.filterFeatures <- function(expr, method) {

  ## filter out genes with constant expression values
  sdGenes <- apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  } 

  expr
}

## maps gene sets content in 'gsets' to 'features', where 'gsets'
## is a 'list' object with character string vectors as elements,
## and 'features' is a character string vector object. it assumes
## features in both input objects follow the same nomenclature,
.mapGeneSetsToFeatures <- function(gsets, features) {

  ## Aaron Lun's suggestion at
  ## https://github.com/rcastelo/GSVA/issues/39#issuecomment-765549620
  gsets2 <- CharacterList(gsets)
  mt <- match(gsets2, features)
  mapdgenesets <- as.list(mt[!is.na(mt)])

  if (length(unlist(mapdgenesets, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  mapdgenesets
}
