# TODO Report where within-study sites have the same coords, dates, 
#      predominant habitat and use intensity

CuratonReport <- function(sites, taxa, bib, diversity) {
    # Takes database extracts and returns a curation as a data.frame.
    old.fragmentation <- c('Part of unfragmented habitat',
                           'Representative part of a fragmented habitat',
                           'Fragment surrounded by other habitats',
                           'Located on the edge of a large continuous habitat',
                           'Overlaps with the intersection of 2 or more continuous habitats')
    old.habitat <- c('Secondary forest',
                     'Secondary non-forest')

    # TODO These checks belong in a EnsureConsistent function?
    stopifnot(all(sites$Source_ID %in% bib$Source_ID))
    stopifnot(all(taxa$Source_ID %in% bib$Source_ID))
    stopifnot(all(diversity$Source_ID %in% bib$Source_ID))
    bib$Source_ID <- as.factor(bib$Source_ID)
    sites$Source_ID <- as.factor(sites$Source_ID)
    taxa$Source_ID <- as.factor(taxa$Source_ID)
    levels(sites$Source_ID) <- levels(bib$Source_ID)
    levels(taxa$Source_ID) <- levels(bib$Source_ID)
    levels(diversity$Source_ID) <- levels(bib$Source_ID)

    curation <- data.frame(
      Invalid_diversity_metric=0!=table(sites[!sites$Diversity_metric_is_valid,'Source_ID']),
      Diversity_metric_requires_curation=0!=table(sites[sites$Diversity_metric_requires_curation,'Source_ID']),
      Missing_sampling_effort=0!=table(sites[sites$Diversity_metric_is_effort_sensitive & (is.na(sites$Sampling_effort) | 0==sites$Sampling_effort),'Source_ID']),
      Invalid_sampling_method=0!=table(sites[!sites$Sampling_method_is_valid,'Source_ID']),
      Sampling_method_requires_curation=0!=table(sites[sites$Sampling_method_requires_curation,'Source_ID']),
      Old_predominant_habitat=0!=table(sites[sites$Predominant_habitat %in% old.habitat,'Source_ID']),
      Old_fragmentation_layout=0!=table(sites[sites$Fragmentation_layout %in% old.fragmentation,'Source_ID']), 
      Old_site_area=0!=table(sites[!is.na(sites$Site_area), 'Source_ID']), 
      Missing_sampling_target=0!=table(sites[''==sites$Sampling_target,'Source_ID']),
      Missing_taxon_status=0!=table(taxa[''==taxa$Status,'Source_ID']), 
      Missing_habitat_area_units=0!=table(MissingHabitatPatchAreaUnits(sites)$Source_ID),
      Bad_blocks=0!=table(BadSiteBlocks(sites)$Source_ID),
      Non_integer=0!=table(NonIntegerValues(diversity)$Source_ID),
      Pseudo_replicated_site_blocks=0!=table(PseudoReplicatedSiteBlocks(sites)$Source_ID),
      Predominant_habitat_cannot_decide=0!=table(sites['Cannot decide'==sites$Predominant_habitat,'Source_ID']), 
      Use_intensity_cannot_decide=0!=table(sites['Cannot decide'==sites$Use_intensity,'Source_ID']), 
      Fragmentation_layout_cannot_decide=0!=table(sites['Cannot decide'==sites$Fragmentation_layout,'Source_ID']), 
      Country_outside_shape=0!=table(sites[sites$Country_distance_metres>0, 'Source_ID']), 
      Eco_region_outside_shape=0!=table(sites[sites$Eco_region_distance_metres>0, 'Source_ID']), 
      Lacking_blocks=0!=table(LackingSiteBlocks(sites)$Source_ID),
      Empty_cells_and_zeroes=0!=table(EmptySiteBySpeciesCells(diversity, TRUE)$Source_ID),
      Empty_cells_without_zeroes=0!=table(EmptySiteBySpeciesCells(diversity, FALSE)$Source_ID),
      Multiple_kingdoms=0!=table(MultipleKingdoms(taxa)$Source_ID),
      Starts_before_2000=0!=table(sites[as.Date(sites$Sample_start_earliest)<as.Date('2000-01-01'), 'Source_ID']),
      Ends_before_2000=0!=table(sites[as.Date(sites$Sample_end_latest)<as.Date('2000-01-01'), 'Source_ID']), 
      Duplicated_taxa=0!=table(DuplicatedTaxa(taxa)$Source_ID), 
      No_gradient=0!=table(NoGradient(sites)$Source_ID),
      All_sites_same_coordinates=0!=table(AllSitesSameCoordinates(sites)$Source_ID), 
      Bad_agricultural=0!=table(BadAgricultural(sites)$Source_ID))

    # Order curation by same Source ID
    curation <- curation[as.character(bib$Source_ID),]
    curation <- cbind(bib[,c('Source_ID',
                             'Curator',
                             'Imported_on', 
                             'Data_format',
                             'Data_use_permission', 
                             'Insightly_permission_email_ID',
                             'Insightly_project', 
                             'Insightly_category',
                             'Insightly_status', 
                             'N_studies',
                             'N_sites',
                             'N_samples',
                             'Sensitive_data')], 
                      curation)

    return (curation)
}

DuplicatedTaxa <- function(taxa) {
    .Log('Computing duplicated taxa\n')
    # Takes a taxa extract and returns a data.frame of studies that contain 
    # taxa with a rank of species or infraspecies that are duplicated.
    taxa <- taxa[,c('Source_ID', 'SS', 'Study_number', 'Taxon_name_entered', 
                    'COL_justification', 'Taxon', 'Rank')]
    taxa <- taxa[taxa$Rank %in% c('Species','Infraspecies'),]
    res <- by(taxa, taxa$SS, function(rows) {
        counts <- table(rows$Taxon)
        duplicated <- names(counts[counts>1])
        return (rows[rows$Taxon %in% duplicated,])
    })

    # Discard non-NULL results with no duplicates
    res <- res[!sapply(res, is.null)]
    res <- res[0!=sapply(res, nrow)]

    res <- do.call('rbind', res)
    rownames(res) <- NULL

    res <- res[order(res$Source_ID, res$Study_number, res$Taxon, res$Taxon_name_entered),]
    return (res)
}

DifferentResolutionsOfSameTaxon <- function(taxa) {
    # Takes a taxa extract and returns a data.frame containing taxon names that 
    # have been resolved to different COL records.
    .Log('Computing different resolutions of the same taxon\n')
    taxa <- taxa[,c('Source_ID', 'Study_number', 'Taxon_name_entered', 
                    'Parsed_name', 'COL_query_name', 'COL_expected_higher_rank', 
                    'COL_expected_higher_taxon')]

    .Log(paste(nrow(taxa), 'rows\n'))

    # Step 1 - consider only non-empty names that appear more than once
    freq <- table(taxa[''!=taxa$Parsed_name,'Parsed_name'])
    consider <- names(freq[freq>1])
    .Log(paste(length(consider), 'non-empty names appear more than once\n'))

    taxa <- droplevels(taxa[taxa$Parsed_name %in% consider,])
    .Log(paste('Considering just these names gives', nrow(taxa), 'rows\n'))

    # Step 2 - count the number of different query names for each parsed name
    res <- tapply(taxa$COL_query_name, taxa$Parsed_name, function(n) {
        return (length(unique(n)))
    })

    res <- res[res>1]
    .Log(paste(length(res), 'names resolved differently\n'))

    # Step 3 - assemble a table showing different resolutions of the same parsed name
    res <- lapply(names(res), function(n) {
        return (taxa[taxa$Parsed_name==n,])
    })

    res <- droplevels(do.call('rbind', res))
    rownames(res) <- NULL
    return (res)
}

EmptySiteBySpeciesCells <- function(diversity, with.zeroes) {
    # Sources with empty cells and zeroes in the site x species matrix
    if(with.zeroes) {
      .Log('Computing empty site by species cells (including zeroes)\n')
    } else {
      .Log('Computing empty site by species cells (excluding zeroes)\n')
    }
    diversity <- diversity[,c('Source_ID', 'Study_number', 'Site_name', 
                              'Taxon_name_entered','Measurement')]
    res <- by(diversity, paste(diversity$Source_ID, diversity$Study_number), 
      function(rows) {
        NTaxa <- length(unique(rows$Taxon_name_entered))
        NSites <- length(unique(rows$Site_name))
        ExpectedNSamples <- NTaxa*NSites
        ActualNSamples <- nrow(rows)
        NEmpty <- ExpectedNSamples-ActualNSamples
        NZeroes <- sum(0==rows$Measurement)
        Sane <- ActualNSamples<=ExpectedNSamples && NEmpty<=ExpectedNSamples
        return (data.frame(Source_ID=rows$Source_ID[1], 
                           Study_number=rows$Study_number[1], 
                           ExpectedNSamples, 
                           ActualNSamples, 
                           NEmpty, 
                           NZeroes, 
                           Sane))
    })

    res <- do.call('rbind', res)

    stopifnot(all(res$Sane))

    if(with.zeroes) res <- res[res$NEmpty>0 & res$NZeroes>0,]
    else            res <- res[res$NEmpty>0 & res$NZeroes==0,]
    rownames(res) <- NULL
    return (res)
}

MultipleKingdoms <- function(taxa) {
    # Takes a taxa extract and returns a data.frame of studies that examine 
    # more than one kingdom
    .Log('Computing multiple kingdoms per study\n')
    taxa <- taxa[,c('Source_ID', 'Study_number', 'Taxon_name_entered', 
                    'Parsed_name', 'COL_query_name', 'Kingdom')]

    # Step 1 - consider only rows with a kingdom
    taxa <- taxa[''!=taxa$Kingdom,]

    # Step 2 - count the number of different kingdoms in each study
    res <- by(taxa, paste(taxa$Source_ID, taxa$Study_number), function(rows) {
        kingdoms <- unique(rows$Kingdom)
        return (data.frame(Source_ID=rows$Source_ID[1], 
                           Study_number=rows$Study_number[1], 
                           N_kingdoms=length(kingdoms),
                           Kingdoms=paste(kingdoms, collapse=',')))
    })

    res <- do.call('rbind', res)

    res <- res[res$N_kingdoms>1,]
    rownames(res) <- NULL
    return (res)
}

LackingSiteBlocks <- function(sites) {
    .Log('Computing studies lacking site blocks\n')
    sites <- sites[,c('Source_ID','SS','Study_number','Block')]
    res <- by(sites, sites$SS, function(rows)
        return (data.frame(Source_ID=rows$Source_ID[1], 
                           Study_number=rows$Study_number[1], 
                           Lacking_blocks=all(''==rows$Block))))
    res <- do.call('rbind', res)
    res <- res[res$Lacking_blocks,]
    return (res)
}

BadSiteBlocks <- function(sites) {
    # Takes a site-level extract and returns a data.frame of studies in which 
    # either some but not all sites have no block or all sites have the same
    # block
    .Log('Computing studies with bad site blocks\n')
    sites <- sites[,c('Source_ID','SS','Study_number','Block')]
    res <- by(sites, sites$SS, function(rows) {
        n <- nchar(as.character(rows$Block))
        ok <- all(0==n) || (all(0!=n) && length(unique(rows$Block))>1)
        return (data.frame(Source_ID=rows$Source_ID[1], 
                           Study_number=rows$Study_number[1],  
                           Blocks_OK=ok))
    })
    res <- do.call('rbind', res)
    res <- res[!res$Blocks_OK,]
    rownames(res) <- NULL
    return (res)
}

PseudoReplicatedSiteBlocks <- function(sites) {
    # Takes a site-level extract and returns a data.frame of studies in which 
    # all sites in a block have the same predominant habitat and use intensity
    .Log('Computing pseudo replicated site blocks\n')
    sites <- sites[,c('Source_ID','SSB','Study_number','Block','Predominant_habitat',
                      'Use_intensity')]
    sites <- sites[''!=sites$Block,]
    res <- by(sites, sites$SSB, 
      function(rows) {
        pr <- 1==nrow(unique(rows[,c('Predominant_habitat','Use_intensity')]))
        return (data.frame(Source_ID=rows$Source_ID[1], 
                           Study_number=rows$Study_number[1],  
                           Block=rows$Block[1],  
                           Pseudo_replicated_blocks=pr))
    })
    res <- do.call('rbind', res)
    res <- res[res$Pseudo_replicated_blocks,]
    rownames(res) <- NULL
    return (res)
}

NoGradient <- function(sites) {
  # Takes a site-level extract and returns a data.frame of studies in which 
  # all sites have the same habitat, use intensity and fragmentation layout.
  .Log('Computing studies with no change in habitat or in use intensity\n')
  res <- by(sites[,c('Source_ID','Study_number','Predominant_habitat',
                     'Use_intensity','Fragmentation_layout')],
            sites$SS, unique)
  res <- res[sapply(res, function(rows) 1==nrow(rows))]
  res <- do.call('rbind', res)
  rownames(res) <- NULL
  return (res)
}

AllSitesSameCoordinates <- function(sites) {
  # Takes a site-level extract and returns a data.frame of studies in which 
  # all sites have the same coordinates.

  # TODO consider precision?
  .Log('Computing studies with all sites at the same coordinates\n')
  res <- by(sites[,c('Source_ID','Study_number','Latitude','Longitude')], 
            sites$SS, unique)
  res <- res[sapply(res, function(rows) 1==nrow(rows))]
  res <- do.call('rbind', res)
  rownames(res) <- NULL
  return (res)
}

MissingHabitatPatchAreaUnits <- function(sites) {
  # Takes a site-level extract and returns a data.frame of studies in which 
  # one or more sites have a habitat patch size and no units.

  .Log('Computing studies with one or more sites with missing habitat patch area units\n')
  res <- by(sites[,c('Source_ID','Study_number','Habitat_patch_area','Habitat_patch_area_unit')], 
    sites$SS, 
    function(rows) {
      bad <- any(!is.na(rows$Habitat_patch_area) & -1!=rows$Habitat_patch_area & ''==rows$Habitat_patch_area_unit)
      return (data.frame(Source_ID=rows$Source_ID[1], 
                         Study_number=rows$Study_number[1],  
                         Missing_habitat_patch_area_unit=bad))
  })
  res <- do.call('rbind', res)
  res <- res[res$Missing_habitat_patch_area_unit,]
  rownames(res) <- NULL
  return (res)
}

BadAgricultural <- function(sites) {
  .Log('Computing studies with invalid agricultural fields\n')

  agro <- sites[,c('Source_ID','SS','Predominant_habitat','Crop','FF1','FF2','FF3','AES','Organic')]

  # Fix NAs
  agro$Crop <- as.character(agro$Crop)
  agro$Crop[is.na(agro$Crop)] <- ''
  agro$Crop <- factor(agro$Crop)
  agro$Organic <- as.character(agro$Organic)
  agro$Organic[is.na(agro$Organic)] <- ''
  agro$Organic['N/A'==agro$Organic] <- ''
  agro$Organic <- factor(agro$Organic)
  agro$AES <- as.character(agro$AES)
  agro$AES[is.na(agro$AES)] <- ''
  agro$AES['N/A'==agro$AES] <- ''
  agro$AES <- factor(agro$AES)

  n <- rowSums(cbind(''!=agro$FF1, ''!=agro$FF2, ''!=agro$FF3, ''!=agro$AES, ''!=agro$Organic))
  crop <- agro[''!=agro$Crop,]
  F <- function(rows, msg) {
    if(nrow(rows))  cbind(unique(rows), msg)
    else            NULL
  }
  res <- rbind(F(agro['Plantation forest'==agro$Predominant_habitat & !substring(agro$FF1,1,1) %in% c('1','3','9'),], 'Incorrect or missing FF1 for Plantation forest'), 
               F(agro['Cropland'==agro$Predominant_habitat & !substring(agro$FF1,1,1) %in% c('1','3','9'),], 'Incorrect or missing FF1 for Cropland'), 
               F(agro['Pasture'==agro$Predominant_habitat & !substring(agro$FF1,1,1) %in% c('2','3','9'),], 'Incorrect or missing FF1 for Pasture'), 
               F(agro[''!=agro$FF1 & !agro$Predominant_habitat %in% c('Plantation forest','Pasture','Cropland'),], 'Incorrect Predominant_habitat for FF1'), 
               F(agro[n>0 & n<5,], 'Some but not all agricultural fields given'), 
               F(crop[!crop$Predominant_habitat %in% c('Cropland','Plantation forest'),], 'Crop given but unexpected Predominant_habitat'), 
               F(crop[''==crop$FF1 | ''==crop$FF2 | ''==crop$FF3,], 'Crop given but not all FF1, FF2 and FF3'), 
               F(crop['1'!=substr(crop$FF3, 1, 1) & !substr(crop$FF3, 1, 3) %in% c('3.1','3.3'),], 'FF3 indicates no crop but Crop was given'), 
               F(agro[!agro$Organic %in% c('','Yes','No','Unknown'),], 'Bad value for Organic'), 
               F(agro[!agro$AES %in% c('','Yes','No','Unknown'),], 'Bad value for AES'))
  res <- res[order(res$Source_ID, res$SS, res$msg),]
  rownames(res) <- NULL
  return (res)
}

NonIntegerValues <- function(diversity) {
  .Log('Computing studies with units that require integer values with one or', 
       'more non-integers values\n')

  require.integers <- c('individuals', 'number of groups', 'number of plots', 
    'number of sites', 'number of species', 'signs', 'species', 
    'times observed')
  diversity <- diversity[,c('Source_ID','SS','Study_number','Measurement','Diversity_metric_unit')]
  diversity <- diversity[diversity$Diversity_metric_unit %in% require.integers,]
  res <- by(diversity, diversity$SS, function(rows) {
      # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
      return (data.frame(Source_ID=rows$Source_ID[1], 
                         Study_number=rows$Study_number[1],  
                         Diversity_metric_unit=rows$Diversity_metric_unit[1], 
                         N_values=nrow(rows), 
                         N_non_integer_values=sum(floor(rows$Measurement)!=rows$Measurement)))
  })
  res <- do.call('rbind', res)
  res <- res[res$N_non_integer_values>0,]
  rownames(res) <- NULL
  return (res)
  res <- data.frame(Source_ID=bad)
}
