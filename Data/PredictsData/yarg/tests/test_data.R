TestTaxonomicCounts <- function() {
  data(taxonomic.groups)
  AssertFalse(any(duplicated(taxonomic.groups$Group)))
  AssertTrue(all(''!=taxonomic.groups$Resolution))
  AssertTrue(all(''!=taxonomic.groups$Higher_group))
  AssertTrue(all(0<taxonomic.groups$N_estimated | is.na(taxonomic.groups$N_estimated)))
  AssertTrue(all(0<taxonomic.groups$N_described | is.na(taxonomic.groups$N_described)))
}
