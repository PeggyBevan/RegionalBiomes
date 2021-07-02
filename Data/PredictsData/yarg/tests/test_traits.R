TestTraitNumeric <- function() {
  res <- TraitDataset('Test', 'Trait 1', 'x', taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=1:2)
  AssertEqual(res$Reference,'Test')
  AssertEqual(res$Trait, 'Trait 1')
  AssertEqual(res$Unit, 'x')
  expected.values <- data.frame(Taxon=c('A','B'), Rank='Species', Kingdom='Animalia', Value=1:2)
  # Factor levels for rank
  expected.values$Rank <- EnsureFactorLevels(expected.values$Rank, TaxonomicRanks())
  AssertEqual(res$Values, expected.values)
}

TestTraitCategorical <- function() {
  res <- TraitDataset('Test', 'Trait 1', 'x', taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=c('CARNIVORE', 'HERBIVORE'))
  AssertEqual(res$Reference,'Test')
  AssertEqual(res$Trait, 'Trait 1')
  AssertEqual(res$Unit, 'x')
  expected.values <- data.frame(Taxon=c('A','B'), Rank='Species', Kingdom='Animalia', Value=c('CARNIVORE', 'HERBIVORE'))
  # Factor levels for rank
  expected.values$Rank <- EnsureFactorLevels(expected.values$Rank, TaxonomicRanks())
  AssertEqual(res$Values, expected.values)
}

TestTraitLogical <- function() {
  res <- TraitDataset('Test', 'Trait 1', 'x', taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=c(TRUE, FALSE))
  AssertEqual(res$Reference,'Test')
  AssertEqual(res$Trait, 'Trait 1')
  AssertEqual(res$Unit, 'x')
  expected.values <- data.frame(Taxon=c('A','B'), Rank='Species', Kingdom='Animalia', Value=c(TRUE, FALSE))
  # Factor levels for rank
  expected.values$Rank <- EnsureFactorLevels(expected.values$Rank, TaxonomicRanks())
  AssertEqual(res$Values, expected.values)
}

TestTraitFailures <- function() {
  # Empty Reference
  AssertRaises(TraitDataset('',     'Trait 1', 'x', taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=1:2))

  # Empty Trait
  AssertRaises(TraitDataset('Test', '',        'x', taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=1:2))

  # Empty Unit
  AssertRaises(TraitDataset('Test', 'Trait 1', '',  taxa=c('A','B'), kingdoms='Animalia', ranks='Species', values=1:2))

  # Duplicated taxa
  AssertRaises(TraitDataset('Test', 'Trait 1', 'x', taxa=c('A','A'), kingdoms='Animalia', ranks='Species', values=1:2))

  # Unknown rank
  AssertRaises(TraitDataset('Test', 'Trait 1', 'x', taxa=c('A','B'), kingdoms='Animalia', ranks='x',       values=1:2))
}

TestAddTraits <- function() {
  input <- read.csv('test_traits_input.csv')
  input$Genus <- sapply(strsplit(as.character(input$Taxon),  ' '), '[', 1)
  traits <- read.csv('test_traits_numeric.csv')
  numeric1 <- with(traits['Numeric 1'==traits$Trait,], TraitDataset('Test','Numeric 1','x',Taxon,'Species',Kingdom,Value))
  numeric2 <- with(traits['Numeric 2'==traits$Trait,], TraitDataset('Test','Numeric 2','x',Taxon,'Species',Kingdom,Value))
  numeric3 <- with(traits['Numeric 3'==traits$Trait,], TraitDataset('Test','Numeric 3','x',Taxon,'Species',Kingdom,Value))

  traits <- read.csv('test_traits_categorical.csv')
  categorical1 <- TraitDataset('Test','Categorical 1','x',traits$Taxon,'Species',traits$Kingdom,traits$Value)

  traits <- read.csv('test_traits_logical.csv')
  logical1 <- TraitDataset('Test','Logical 1','x',traits$Taxon,'Species',traits$Kingdom,traits$Value)

  res <- input
  res <- AddTraits(res, numeric1, numeric2, numeric3, categorical1, logical1)
  res <- AddGenusAveragedTraits(res, numeric1, numeric2, numeric3)

  # Can't genus-averaged categorical or logical traits
  AssertRaises(AddGenusAveragedTraits(res, categorical1))
  AssertRaises(AddGenusAveragedTraits(res, logical1))

  expected <- read.csv('test_traits_expected.csv', check.names=FALSE)

  res$Categorical_1_x <- as.character(res$Categorical_1_x)
  expected$Categorical_1_x <- as.character(expected$Categorical_1_x)
  AssertEqual(expected[,intersect(colnames(res), colnames(expected))], 
              res[,intersect(colnames(res), colnames(expected))])

  # Existing trait values should not be altered
  res <- AddTraits(res, numeric1, numeric2, numeric3, categorical1)
  res <- AddGenusAveragedTraits(res, numeric1, numeric2, numeric3)
  AssertEqual(expected[,intersect(colnames(res), colnames(expected))], 
              res[,intersect(colnames(res), colnames(expected))])

  cat('TODO test InferTraitValues\n')
}
