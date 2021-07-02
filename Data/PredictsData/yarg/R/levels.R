# Factor levels
EnsureFactorLevels <- function(x, levels, allow.empty=FALSE, ...) {
  if(allow.empty && !'' %in% levels) {
    levels <- c('', levels)
  }

  unrecognised <- setdiff(x, levels)
  if(length(unrecognised)>0) {
    stop('Unrecognised levels ', paste(unrecognised, collapse=','))
  } else {
    return (factor(x, levels, ...))
  }
}

UseIntensities <- function() {
  return (c("Minimal use", 
            "Light use",
            "Intense use",
            "Cannot decide"))
}

Habitats <- function() {
  return (c("Primary forest",
            "Primary non-forest",
            "Young secondary vegetation",
            "Intermediate secondary vegetation",
            "Mature secondary vegetation",
            "Secondary vegetation (indeterminate age)",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban",
            "Cannot decide"))
}

SimplifyHabitat <- function(v) {
  # Combines two primary levels in v
  simplified <- c("Primary vegetation",
                  "Young secondary vegetation",
                  "Intermediate secondary vegetation",
                  "Mature secondary vegetation",
                  "Secondary vegetation (indeterminate age)",
                  "Plantation forest",
                  "Pasture",
                  "Cropland",
                  "Urban",
                  "Cannot decide")
  v <- as.character(v)
  v['Primary forest'==v] <- 'Primary vegetation'
  v['Primary non-forest'==v] <- 'Primary vegetation'
  v <- EnsureFactorLevels(v, simplified)
  return (v)
}

FragmentationLayouts <- function() {
  return (c('Well within unfragmented habitat',
            'Within unfragmented habitat but at or near its edge',
            'Within remnant patch (perhaps at its edge) that is surrounded by other habitats',
            'Representative part of a fragmented landscape',
            'Part of the matrix surrounding remnant patches',
            'Cannot decide', 
            'Not evaluated'))
}

Biomes <- function() {
  data(biomes)
  return (as.character(biomes$Biome))
}

Ecoregions <- function() {
  data(ecoregions)
  return (sort(as.character(ecoregions$Ecoregion)))
}

Countries <- function() {
  data(countries)
  return (sort(as.character(countries$Country)))
}

UNRegions <- function() {
  data(countries)
  return (sort(as.character(unique(countries$UN_region))))
}

UNSubRegions <- function() {
  data(countries)
  return (sort(as.character(unique(countries$UN_subregion))))
}

Hotspots <- function() {
  data(hotspots)
  return (sort(as.character(hotspots$Hotspot)))
}

WildernessAreas <- function() {
  data(wilderness.areas)
  return (sort(as.character(wilderness.areas$Wilderness_area)))
}

Realms <- function() {
  data(realms)
  return (sort(as.character(realms$Realm)))
}

TaxonomicRanks <- function() {
  return (c('Infraspecies','Species','Genus','Family','Order','Class','Phylum',
            'Kingdom'))
}

DiversityMetricTypes <- function() {
  return (c('Abundance','Occurrence','Species richness'))
}

FF1 <- function() {
  return (c('1. Plant crops', '2. Livestock', '3. Mixed farming', '9. Unknown'))
}

Organic <- function() {
  return (c('Yes','No','Unknown'))
}

AES <- function() {
  return (c('Yes','No','Unknown'))
}

IUCNStatuses <- function() {
  # A mix of 1994 and 2001 statuses
  return (c('EX','EW','CR','EN','VU','NT','LC','LR/cd','LR/nt','LR/lc','DD','NE'))
}

SimplifyIUCN <- function(v) {
  # Encode IUCN Red List status as number between 0 and 5, following Purvis et
  # al 2000 Proc R Soc B
  mapping <- c(EX=5,
               EW=5,
               CR=4,
               EN=3,
               VU=2,
               'LR/cd'=2,
               LC=1,
               'LR/nt'=1,
               'LR/lc'=1,
               NT=0)
  stopifnot(all(is.na(v) | 'DD' == v | !is.na(match(v, names(mapping)))))
  return (mapping[v])
}
