TestSiteMetrics <- function() {
  # Numeric
  trait1 <- TraitDataset(reference='Test', trait='M', unit='kg', 
                         taxa=paste('Sp', 1:6), kingdoms='Animalia', 
                         ranks='Species', values=1:6)

  # Logical
  trait2 <- TraitDataset(reference='Test', trait='Specialist', unit='Logical', 
                         taxa=paste('Sp', 1:5), kingdoms='Animalia', 
                         ranks='Species', values=c(FALSE,TRUE,TRUE,TRUE,TRUE))

  # TODO Categorical

  diversity <- read.csv('test_site_metrics_input.csv')
  diversity <- AddTraits(diversity, trait1, trait2)

  # Add SS and SSS
  diversity <- AddSummaryColumns(diversity)

  # Preprocessing
  diversity <- DropInvalidMetricsAndMethods(diversity)
  diversity <- CorrectSamplingEffort(diversity)

  actual <- SiteMetrics(diversity, traits=c('M_kg', 'Specialist_Logical'))

  cat('HACK: TestSiteMetrics ignoring Richness_rarefied.\n')
  expected <- read.csv('test_site_metrics_expected.csv')
  actual <- actual[,colnames(expected)]
  AssertEqual(expected, actual)

  cat('HACK: TestSiteMetrics setting Taxon_name_entered for non-unique sites.\n')
  diversity$Taxon_name_entered <- diversity$Taxon
  actual <- SiteMetrics(diversity, traits=c('M_kg', 'Specialist_Logical'), sites.are.unique=FALSE)
  actual <- actual[,colnames(expected)]
  AssertEqual(expected, actual)

  # Can't compute site metrics of categorical traits
  categorical <- TraitDataset(reference='Test', trait='Diet', unit='Category', 
                              taxa=paste('Sp', 1:3), kingdoms='Animalia', 
                              ranks='Species', 
                              values=c('Herbivore', 'Omnivore', 'Carnivore'))
  diversity <- AddTraits(diversity, categorical)
  AssertRaises(SiteMetrics(diversity, traits='Diet_Category'))
}
