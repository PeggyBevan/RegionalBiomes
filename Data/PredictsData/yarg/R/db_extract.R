SetFactorLevels <- function(extract) {
  # Set factor levels
  F <- function(col, fn, allow.empty=FALSE) {
    if(col %in% colnames(extract)) {
      extract[,col] <<- EnsureFactorLevels(extract[,col], fn(), allow.empty)
    }
  }

  F('Rank', TaxonomicRanks, TRUE) # Uncurated names might not have a Rank
  F('Predominant_habitat', Habitats)
  F('Use_intensity', UseIntensities)
  F('Biome', Biomes)
  F('Ecoregion', Ecoregions)
  F('Country', Countries)
  F('UN_region', UNRegions)
  F('UN_subregion', UNSubRegions)
  F('Hotspot', Hotspots, TRUE)    # Not all sites will have a hotspot
  F('Wilderness_area', WildernessAreas, TRUE)    # Not all sites will have a wilderness area
  F('Realm', Realms)
  F('Diversity_metric_type', DiversityMetricTypes)
  F('Fragmentation_layout', FragmentationLayouts)

  if('FF1' %in% colnames(extract)) {
    # Fix old definitions if required
    from <- c('1. ARABLE',     '2. LIVESTOCK','3. MIXED FARMING')
    to <-   c('1. Plant crops','2. Livestock','3. Mixed farming')
    if(any(from %in% levels(extract$FF1))) {
      .Log('Correcting old FF1 levels\n')
      extract$FF1 <- factor(extract$FF1, union(levels(extract$FF1), FF1()))
      for(i in 1:length(from)) {
        extract$FF1[extract$FF1==from[i]] <- to[i]
      }
    }
    F('FF1', FF1, TRUE)    # Not all sites will have an FF classification
  }

  F('AES', AES, TRUE)    # Not all sites will have an AES classification
  F('Organic', Organic, TRUE)    # Not all sites will have an organic classification

  return (extract)
}

AddSummaryColumns <- function(extract) {
  # Adds columns SS, SSS, SSB and SSBS to extract, if not present and if 
  # possible.

  # TODO This depends upon diversity being ordered by Source_ID, Study_number 
  # and Site_number - enforce or check this ordering?

  F <- function(extract, col, from) {
    # A helper function that adds a new column 'col' to extract. Values in col 
    # will be values in columns given in 'from' pasted together.
    if(!col %in% colnames(extract) && all(from %in% colnames(extract))) {
      .Log('Adding [', col, ']\n', sep='')
      # levels is included in the call to factor() to prevent factor() from 
      # ordering levels alphabetically, which is silly because we want to 
      # retain numerical ordering of study numbers, site numbers and block 
      # numbers. 
      vals <- do.call('paste', extract[,from])
      extract[,col] <- factor(vals, levels=unique(vals))
    }
    return (extract)
  }

  extract <- F(extract, 'SS',   c('Source_ID','Study_number'))
  extract <- F(extract, 'SSS',  c('Source_ID','Study_number', 'Site_number'))
  extract <- F(extract, 'SSB',  c('Source_ID','Study_number', 'Block'))
  extract <- F(extract, 'SSBS', c('Source_ID','Study_number', 'Block', 'Site_number'))

  return (extract)
}

AddStudyCommonTaxon <- function(extract) {
  # Returns extract with two columns added: Study_common_taxon and 
  # Rank_of_study_common_taxon

  ranks <- head(rev(TaxonomicRanks()), -1)
  cols <- c(ranks, c('Higher_taxon','Source_ID','Study_name','SS'))

  add <- c('Study_common_taxon','Rank_of_study_common_taxon')

  if(!any(add %in% colnames(extract)) && all(cols %in% colnames(extract))) {
    .Log('Computing study common taxon for', length(unique(extract$SS)), 
         'studies in', nrow(extract), 'rows\n')

    ct <- by(extract[,cols], 
      extract$SS, 
      function(rows) {
        taxon <- rank <- ''
        for(r in ranks) {
          consider.rank <- FALSE
          if('Species'==r && !all(''==rows[,'Species'])) {
            common <- unique(paste(rows[,'Genus'],rows[,'Species']))
            consider.rank <- TRUE
          } else if('Species'!=r) {
            common <- unique(rows[,r])
            consider.rank <- TRUE
          }

          if(consider.rank && 1==length(common) && ''!=common) {
            taxon <- common
            rank <- r
          } else {
            break;
          }
        }

      return (data.frame(Source_ID=rows$Source_ID[1], 
                         Study_name=rows$Study_name[1],
                         SS=rows$SS[1], 
                         Study_common_taxon=taxon, 
                         Rank_of_study_common_taxon=rank))
    })

    ct <- do.call('rbind', ct)
    .Log(paste('Adding [', add, ']', sep=''), sep='\n')
    extract[,add] <- ct[match(extract$SS, ct$SS),add]

    # Set factor levels
    extract$Study_common_taxon <- droplevels(extract$Study_common_taxon)
    extract$Rank_of_study_common_taxon <-
                      EnsureFactorLevels(extract$Rank_of_study_common_taxon,
                                         TaxonomicRanks(), allow.empty=TRUE)
  }

  return (extract)
}

AddBestGuessBinomial <- function(extract) {
  # Returns extract with Best_guess_binomial added:
  # Taxon name if Rank is Species
  # First two words of Taxon name if Rank is Infraspecies
  # Parsed_name if Rank is neither Species nor Infraspecies and Parsed_name 
  # contains two words.
  required <- c('Taxon', 'Rank', 'Parsed_name')
  if(!'Best_guess_binomial' %in% colnames(extract) && 
     all(required %in% colnames(extract))) {
    .Log('Adding [Best_guess_binomial]\n')

    extract$Best_guess_binomial <- ''

    # Use species names if present
    i <- 'Species'==extract$Rank
    extract$Best_guess_binomial[i] <- as.character(extract$Taxon[i])

    # First two words of subspecies
    empty <- ''==extract$Best_guess_binomial
    if(any(empty)) {
      i <- 'Infraspecies'==extract$Rank
      re <- '^(\\w+) (\\w+) (\\w+)$'
      extract$Best_guess_binomial[i] <- gsub(re, '\\1 \\2', extract$Taxon[i], perl=TRUE)
    }

    # First two words if Parsed_name looks like a binomial or an infraspecies
    empty <- ''==extract$Best_guess_binomial
    if(any(empty)) {
      re <- '^(\\w+) (\\w+)( \\w+)?$'
      i <- grep(re, extract$Parsed_name[empty], perl=TRUE)
      v <- gsub(re, '\\1 \\2', extract$Parsed_name[empty][i], perl=TRUE)
      v <- Capitalize(v)
      extract$Best_guess_binomial[empty][i] <- v
    }

    extract$Best_guess_binomial <- factor(extract$Best_guess_binomial)
  }

  return (extract)
}

AddingSampleMidpoint <- function(extract) {
    # Adds the 'Sample_midpoint' column if extract contains
    # Sample_start_earliest and Sample_end_latest
    if('Date' == class(extract$Sample_start_earliest) &&
       'Date' == class(extract$Sample_end_latest)) {
      extract$Sample_midpoint <- extract$Sample_start_earliest +
                    (extract$Sample_end_latest - extract$Sample_start_earliest) / 2
    }
    return (extract)
}

ReadDBExtract <- function(file, nrows=-1, parse.dates=FALSE, ...) {
    # Reads a CSV extract from the database.
    # Sets factor levels as appropriate.

    # TODO ReadDBExtract to validate DB extracts

    # TODO Much code depends upon diversity being ordered by Source_ID, 
    # Study_number, Site_number, Taxon_number - enforce or check this ordering?

    # Classes for columns in the full diversity-level extracts.
    col.classes <- c(Source_ID='factor',
                     Data_use_permission='factor',
                     Study_number='integer',
                     Study_name='factor',
                     Extent='factor',
                     Location='factor',
                     Region='factor',
                     Sampling_target='factor',
                     Diversity_metric='factor',
                     Diversity_metric_unit='factor',
                     Diversity_metric_type='factor',
                     Diversity_metric_is_valid='logical',
                     Diversity_metric_requires_curation='logical',
                     Diversity_metric_is_effort_sensitive='logical',
                     Sampling_method='factor',
                     Sampling_effort_unit='factor',
                     Sampling_method_is_valid='logical',
                     Sampling_method_requires_curation='logical',
                     Site_number='integer',
                     Site_name='factor',
                     Block='factor',
                     Sample_date_resolution='factor',
                     Country_distance_metres='numeric',
                     Coordinates_method='factor',
                     Coordinates_precision='numeric',
                     Coordinates_precision_unit='factor',
                     Site_area='numeric', 
                     Site_area_unit='factor',
                     Sample_effort='numeric', 
                     Habitat_as_described='factor',
                     Predominant_habitat='factor',
                     Source_for_predominant_habitat='factor',
                     Use_intensity='factor',
                     Fragmentation_layout='factor',
                     Km_to_nearest_edge_of_habitat='numeric',
                     Years_since_fragmentation_or_conversion='numeric',
                     Transect_details='factor',
                     FF1='factor',
                     FF2='factor',
                     FF3='factor',
                     AES='factor',
                     Organic='factor',
                     Crop='factor',
                     Eco_region_distance_metres='numeric',
                     Longitude='numeric',
                     Latitude='numeric',
                     Country='factor',
                     UN_subregion='factor',
                     UN_region='factor',
                     Ecoregion='factor',
                     Biome='factor',
                     Realm='factor',
                     Hotspot='factor',
                     Wilderness_area='factor',
                     Taxon_number='integer',
                     Taxon_name_entered='factor',
                     Resolution_entered='factor',
                     Indication='factor',
                     Status='factor',
                     COL_ID='integer',
                     Taxon='factor',
                     Name_status='factor',
                     Rank='factor',
                     Kingdom='factor',
                     Phylum='factor',
                     Class='factor',
                     Order='factor',
                     Family='factor',
                     Genus='factor',
                     Species='factor',
                     Higher_taxon='factor',
                     Measurement='numeric')

    if(parse.dates) {
        # Parsing dates is slow so is optional

        date.cols <- c('Compiled_on',
                       'Insightly_dont_publish_before',
                       'Sample_start_earliest',
                       'Sample_end_latest')
        cols <- rep('DBExtractDate', length(date.cols))
        names(cols) <- date.cols
        col.classes <- c(col.classes, cols)

        datetime.cols <- c('Added_on',
                           'COL_query_date',
                           'Imported_on')
        cols <- rep('DBExtractDateTime', length(datetime.cols))
        names(cols) <- datetime.cols
        col.classes <- c(col.classes, cols)
    }

    # Read the column-names in the csv file
    cols <- colnames(read.csv(file=file, nrows=1, ...))

    extract <- read.csv(file=file, colClasses=col.classes[cols], nrows=nrows, ...)
    extract <- SetFactorLevels(extract) 
    extract <- AddSummaryColumns(extract)
    extract <- AddStudyCommonTaxon(extract)
    extract <- AddBestGuessBinomial(extract)
    extract <- AddingSampleMidpoint(extract)

    return (extract)
}
