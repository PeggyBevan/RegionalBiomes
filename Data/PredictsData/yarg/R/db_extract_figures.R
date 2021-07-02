# Visualise database extracts
BlendColour <- function(a, b, weight=0.5) {
    stopifnot(0 < weight && weight < 1)
    a <- col2rgb(a)
    b <- col2rgb(b)
    lower <- pmin(a, b)
    upper <- pmax(a, b)
    colour <- lower + weight * (upper - lower)
    to.rgb <- function(c) {
        return(rgb(c[1], c[2], c[3], maxColorValue=0xff))
    }
    return(apply(colour, FUN=to.rgb, MARGIN=2))
}

TaxonomicGroupStrongColours <- function() {
  # My original
  return (c(Arthropoda="#FF00FF",      # purple
            Chordata="#FF0000",        # red
            Fungi="#007EFF",           # blue
            Other="grey",
            Plantae="#00FF00"          # green
           ))
}

TaxonomicGroupSoftColours <- function() {
  # Color brewer 2 with my tinkering
  return (c(Arthropoda="#EB79EB", # B683BE
            Chordata="#e41a1c", 
            Fungi="#377eb8", 
            Other="darkgrey",
            Plantae="#4daf4a"))
}

LandUseStrongColours <- function() {
  return (c('Primary vegetation'="#66A61E",
            'Secondary vegetation'="#1B9E77",
            'Wood plantation'="#7570B3",
            Pasture="#E6AB02",
            Cropland="#D95F02",
            Urban="#E7298A"))
}

CountSpecies <- function(taxa) {
  # The number of species
  # Taxa resolved to species level count once per taxon name
  # Higher ranks count once per name per study

  # Take unique combinations of Source_ID, Study_name and Taxon_name_entered 
  # to allow taxa to be a data.frame of diversity measurements.
  taxa <- taxa[!duplicated(taxa[,c('Source_ID','Study_name','Taxon_name_entered')]),]
  return (length(unique(taxa[taxa$Rank %in% c('Species', 'Infraspecies'),'Taxon'])) + 
          nrow(taxa[!taxa$Rank %in% c('Species', 'Infraspecies'),]))
}

CountSpeciesByGroup <- function(taxa, group) {
  # A named vector of the number of species in each group
  # Take just the rows that we are interested in
  taxa <- taxa[,c('Source_ID','Study_name','Taxon_name_entered','Taxon','Rank',group)]

  # Take unique combinations of Source_ID, Study_name and Taxon_name_entered 
  # to allow taxa to be a data.frame of diversity measurements.
  taxa <- taxa[!duplicated(taxa[,c('Source_ID','Study_name','Taxon_name_entered')]),]

  # Taxa resolved to species level count once per taxon name
  species.level <- taxa$Rank %in% c('Species', 'Infraspecies')
  a <- tapply(taxa[species.level,'Taxon'], droplevels(taxa[species.level,group]), 
              function(taxa) length(unique(taxa)))

  # Higher ranks count once per name per study
  b <- table(taxa[!species.level,group], useNA='ifany')

  names(a)[''==names(a) | is.na(names(a))] <- '<unnamed>'
  names(b)[''==names(b) | is.na(names(b))] <- '<unnamed>'

  # names(b) computed from factor levels
  group.names <- names(b)
  res <- rep(0, length(group.names))
  names(res) <- group.names
  res[names(a)] <- a
  res[names(b)] <- res[names(b)] + b

  stopifnot(sum(res)==CountSpecies(taxa))

  return (res)
}

CountSpeciesByHigherGroup <- function(taxa) {
  # A data.frame of the number of species in each higher taxon together 
  # with estimates of global described N species, global number of 
  # species and the name of the higher group.

  # The counts in the database
  N <- CountSpeciesByGroup(taxa, 'Higher_taxon')

  # Alters counts in N by consolidating and/or renaming names in PREDICTS to
  # match names presented by Chapman 2009
  F <- function(N, predicts, chapman) {
      stopifnot(1==length(chapman))
      stopifnot(0<length(predicts))
      already.present <- intersect(chapman, names(N))
      if(length(already.present)) {
          stop(paste('Names already present: ',
                     paste(already.present, collapse=', ')))
      } else {
          N[chapman] <- sum(N[predicts], na.rm=TRUE)
          N <- N[setdiff(names(N), predicts)]
          return (N)
      }
  }

  # Plants
  N <- F(N, c('Magnoliopsida','Liliopsida'), 'Magnoliophyta')
  N <- F(N, c('Cycadopsida', 'Gnetopsida', 'Pinopsida'), 'Gymnosperms')
  N <- F(N, c('Equisetopsida','Lycopodiopsida','Marattiopsida',
              'Polypodiopsida','Psilotopsida'), 'Ferns and allies')

  # Fungoid protists - from Chapman 2009 p 62
  N <- F(N, c('Mycetozoa', 'Myxomycota', 'Plasmodiophoromycota'), 'Fungoid protists')

  # Psocodea "...was formerly considered a superorder, but is now generally
  # considered by entomologists as an order." and are "believed to have evolved
  # from within the former order "Psocoptera", which contained the bark lice and
  # book lice."
  # http://en.wikipedia.org/wiki/Psocodea - 2015-01-14
  N <- F(N, 'Psocodea', 'Psocoptera')

  # Malacostraca and Maxillopoda are the only classes in Crustacea that have 
  # terrestrial species (I think)
  N <- F(N, c('Malacostraca', 'Maxillopoda'), 'Crustacea')

  N <- F(N, 'Zygentoma', 'Thysanura')

  # Sort by increasing counts the alphabetically
  N <- N[order(-N, names(N), decreasing=TRUE)]

  data(taxonomic.groups)
  rownames(taxonomic.groups) <- as.character(taxonomic.groups$Group)

  # Assemble N counts, estimate of numbers of described species and estimate of the total numbers of species
  counts <- cbind(N_represented=N, taxonomic.groups[names(N),])

  # Missing names are included in the returned data.frame but will not have
  # values from taxonomic.groups
  missing <- setdiff(names(N), c('<unnamed>', rownames(taxonomic.groups)))
  if(length(missing)) {
    cat('Unrecognised taxonomic groups', paste(sort(missing), collapse=', '), '\n')
  }

  return (counts)
}

PlotTaxonRanks <- function(taxa, 
                           mar=c(4, 3.2, 1, 0.2), # bottom, left, top, right
                           mgp=c(2, 0.6, 0),      # axis title, axis labels and axis line
                           cummulative=TRUE, label.offset=-1,
                           srt=35, cex.axis=par('cex.axis'), ...) {
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))
  ranks <- CountSpeciesByGroup(taxa, 'Rank')
  ranks <- 100*ranks/sum(ranks)
  if(cummulative) ranks <- cumsum(ranks)
  mids <- barplot(ranks, xaxt='n', ylab='% species', ...)
  # Rotate x axis labels
  # http://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-I-create-rotated-axis-labels_003f
  text(x=mids+(mids[2]-mids[1])/4, y=par("usr")[3]+label.offset, 
       labels=names(ranks), xpd=TRUE, srt=srt, cex=cex.axis, adj=1, ...)
  invisible (mids)
}

.BarText <- function(x, y, labels, right.threshold, left.offset, right.offset, 
                     ...) {
  # You might be tempted to use the pos argument of text() - this has the 
  # side-effect of shifting the text vertically
  right <- x<right.threshold & !is.na(x)
  if(any(right)) {
    text(x[right], y[right], labels[right], adj=c(right.offset,NA), ...)
  }

  left <- !right
  if(any(left)) {
    text(x[left], y[left], labels[left], adj=c(left.offset,NA), ...)
  }
}

HorizontalBarplot <- function(values, 
                              labels=format(values, big.mark=',', big.interval=3), 
                              right.threshold=max(values, na.rm=TRUE)/30, 
                              cex.axis=par("cex.axis"),
                              cex.names=par("cex.axis"), las=1, 
                              show.labels=TRUE, left.offset=1,
                              right.offset=0, ...) {
  # TODO format thousands
  # Do not show ugly single-pixel width lines or text for zeroes. 
  values[0==values] <- NA

  y <- barplot(values, horiz=TRUE, cex.axis=cex.axis, cex.names=cex.names, 
               las=las, ...)
  if(show.labels) {
    .BarText(values, y, labels, right.threshold, left.offset, right.offset, ...)
  }
  invisible(y)
}

PlotTaxonCountsBar <- function(counts, 
                               mar=c(3, 4, 0.2, 0.2), # bottom, left, top, right
                               mgp=c(1.5, 0.6, 0),    # axis title, axis labels and axis line
                               col=TaxonomicGroupSoftColours()[counts$Higher_group], 
                               border=col, yaxs='i', ...) {
  # yaxs defaults to 'i' (internal) to prevent the gap of 4% of the 
  # data's range from being added between the data and the axis. See help page 
  # for par.
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))

  if(any(is.na(border)) && !all(is.na(border))) {
      border[is.na(border)] <- 'lightgrey'
  }

  bar.height <- 0.8
  n <- nrow(counts)
  log10counts <- log10(counts[,c('N_represented','N_described','N_estimated')])

  plot(NA, NA, xlim=c(0, max(log10counts, na.rm=TRUE)), 
       ylim=c(1-bar.height/2,n+bar.height/2), frame.plot=FALSE, 
       yaxt='n', ylab='', xlab=~log[10](N), yaxs=yaxs, ...)
  axis(side=2, at=1:n, labels=rownames(log10counts), tick=FALSE, las=1, pos=0, ...)
  segments(0, 1:n, pmax(log10counts[,2], log10counts[,3], na.rm=TRUE), 1:n, ...)
  points(log10counts[,'N_estimated'], 1:n, pch=19, ...)    # filled circles 
  points(log10counts[,'N_described'], 1:n, pch=4, ...)     # crosses
  rect(0, (1:n)-bar.height/2, log10counts[,'N_represented'], (1:n)+bar.height/2, 
       col=col, border=border, ...)
  # Use cex.axis, if set
  counts.text <- paste(FormatThousands(counts[,'N_represented']), ' (', 
                       sprintf('%.1f%%', 100*counts[,'N_represented']/counts[,'N_described']), ')', 
                       sep='')
  counts.text[is.na(counts[,'N_described'])] <- counts[is.na(counts[,'N_described']),'N_represented']
  dots <- list(...)
  # Do not show text to the right of the bars - it looks ugly on top of the
  # horizontal lines
  counts.text[counts[,'N_represented']<5] <- ''
  .BarText(x=log10counts[,'N_represented'], y=1:n, labels=counts.text,
           right.threshold=0, left.offset=1.1, right.offset=0,
           cex=ifelse('cex.axis' %in% names(dots), dots[['cex.axis']], par('cex.axis')))
}

PlotTaxonCounts <- function(counts, 
                            mar=c(3, 3.2, 0.5, 0.2), # bottom, left, top, right
                            mgp=c(2, 0.6, 0),        # axis title, axis labels and axis line
                            tck=NA,
                            ps= NA,
                            bg=TaxonomicGroupStrongColours()[counts$Higher_group],
                            col=1,
                            pch=21, 
                            label.right=NULL,      # Names of groups for which labels will be shown to the right of the point
                            label.above=NULL,      # Names of groups for which labels will be shown above the point
                            label.below=NULL,      # Names of groups for which labels will be shown below the point
                            labels.cex=0.7, 
                            cex=1.2, cex.legend=1.2,
                            pt.cex=1.2,
                            clever=FALSE,
                            spaces=1,
                            ...) {
  if (is.na(ps)) ps <- par("ps")
  
  op <- par(mar=mar, mgp=mgp,tck=tck,ps=ps)
  on.exit(par(op))

  counts <- counts[!is.na(counts$N_described),]
  x <- log10(counts[,c('N_described')])
  y <- log10(counts[,c('N_represented')])
  plot(x, y, xlab=~Log[10](estimated~number~of~described~species),
       ylab=~Log[10](number~of~species~'in'~database), type='n',
       cex=cex, ...)

  abline(a=log10(0.001), b=1, col='black', lty=3, ...)  #  0.1% line
  abline(a=log10(0.01),  b=1, col='black', ...)         #  1% line
  abline(a=log10(0.10),  b=1, col='black', lty=2, ...)  # 10% line

  if(clever) {
    require(wordcloud)

    # Label positions
    nc <- wordlayout(c(x,x), c(y,y), c(rep('X', nrow(counts)), rownames(counts)),
                     xlim=range(x), ylim=range(y), cex=labels.cex)
    nc <- nc[(1+(nrow(nc)/2)):nrow(nc),]

    # Lines
    segments(x, y, nc[,1]+nc[,3]/2, nc[,2]+nc[,4]/2, col='lightgrey', ...)

    # Points
    points(x, y, col=col, bg=bg, pch=pch, cex=pt.cex, ...)

    # Labels
    text(nc[,1] + .5*nc[,3],nc[,2]+.5*nc[,4], rownames(counts), cex=labels.cex, 
         ...)
  } else {
    points(x, y, col=col, bg=bg, pch=pch, cex=pt.cex, ...)

    # Place labels to the left or right
    label.margin <- paste(rep(' ', spaces), collapse='')

    left <- !rownames(counts) %in% c(label.right, label.above, label.below)
    labels.locations <- list(
      list(show=left, adj=c(1, 0.45), prefx='', suffix=label.margin),
      list(show=rownames(counts) %in% label.right, adj=c(0, 0.45),
           prefix=label.margin, suffix=''),
      list(show=rownames(counts) %in% label.above, adj=c(0.5, -1.1),
           prefix='', suffix=''),
      list(show=rownames(counts) %in% label.below, adj=c(0.5, 2),
           prefix='', suffix='')
    )

    for(row in 1:length(labels.locations)) {
      show <- labels.locations[[row]]$show
      if(any(show)) {
        labels <- paste0(labels.locations[[row]]$prefix,
                         rownames(counts)[show],
                         labels.locations[[row]]$suffix)
        text(x[show], y[show], labels=labels, adj=labels.locations[[row]]$adj,
             cex=labels.cex, ...)
      }
    }
  }

  legend('topleft', legend=c('  0.1%', '  1.0%', '10.0%'), lty=c(3, 1, 2),
    title='Representation', cex=cex.legend)
}

PlotLocationsMap <- function(locations, mar=rep(0,4), ...) {
  require(maps)
  op <- par(mar=mar)
  on.exit(par(op))
  map('world', col='darkgrey', mar=rep(0, 4), lwd=0.5, ...)
  points(locations$Longitude, locations$Latitude, pch=19, ...)
}

PlotLocations <- function(locations, biomes.grid, mar=rep(0,4), cex=1,
                          points.cex, radii, legend.cex=cex, pch=19,
                          bg='#99B3CC', blend.weight=0.5, col='black',
                          box.lwd=2, ...) {
  # locations - a data.frame with columns Latitude and Longitude
  # Either circles or points
  stopifnot(1==(missing(points.cex) + missing(radii)))
  op <- par(mar=mar)
  on.exit(par(op))
  require(rgdal)

  b <- biomes.grid

  # Don't plot empty cells
  clear <- 0==b$band1 & 0==b$band2 & 0==b$band3
  b$band1[clear] <- b$band2[clear] <- b$band3[clear] <- NA

  if(!is.null(blend.weight)) {
    # Blend none-white bands towards white
    stopifnot(0 < blend.weight && blend.weight < 1)
    white <- !clear & 0xff==b$band1 & 0xff==b$band2 & 0xff==b$band3

    BlendValue <- function(a, b)  {
      lower <- pmin(a, b)
      upper <- pmax(a, b)
      return (lower + blend.weight * (upper - lower))
    }

    b$band1[!white] <- BlendValue(b$band1[!white], 0xff)
    b$band2[!white] <- BlendValue(b$band2[!white], 0xff)
    b$band3[!white] <- BlendValue(b$band3[!white], 0xff)
  }

  image(b, bg=bg, ylim=c(-120, 85), red="band1", green="band2", blue="band3")

  if(!missing(points.cex)) {
    points(Latitude~Longitude, data=locations, cex=points.cex,
           bg=col, col=col, pch=21, ...)
  } else {
    symbols(locations$Longitude, locations$Latitude, circles=radii, bg=col,
            fg=NA, add=TRUE, inches=FALSE)
  }

  data(biomes)
  biomes <- biomes[order(biomes$ID),]
  biomes <- biomes[!biomes$Biome %in% c('Inland Water', 'Rock & Ice'),]
  legend(-190, -62, legend=paste(biomes$ID, biomes$Biome), 
         fill=as.character(biomes$Weak_colour), border='black', 
         text.col='black', horiz=FALSE, x.intersp=0.5,
         text.width=175, ncol=2, xpd=TRUE, cex=legend.cex, 
         bg=bg, bty='o', box.col='black', box.lwd=box.lwd, ...)
}

PlotNPerCountryChloropleth <- function(N, world.borders, mar=rep(0,4), ...) {
  # Chloropleth showing N per country
  # N should be a vector of counts. Names of N should be country names.
  op <- par(mar=mar)
  on.exit(par(op))

  require(rgdal)
  stopifnot(all(names(N) %in% world.borders$NAME))

  colours <- c("#D4B9DA", "#C894BB", "#BC6F9D", "#B04A7F", "#A42461")
  colours <- rev(heat.colors(8))
  col <- rep('grey', length(world.borders$NAME))
  names(col) <- as.character(world.borders$NAME)
  breaks <- seq(1,5, length.out=1+length(colours))
  col[names(N)] <- colours[cut(log10(N), breaks=breaks, labels=FALSE)]
  par(mar=rep(0,4))
  plot(world.borders, border="white", lwd=1, col=col)
  legend(-160, -95, fill=colours, ncol=length(colours)/2, bg='white', 
         legend=paste(head(breaks,-1), '<=log10(n)<', head(breaks,-1)+diff(breaks)[1], sep=''))
}

.RegionsTable <- function(regions, col, sites, diversity) {
  # regions - a data.frame
  # col - the name of a column in regions, sites and diversity
  res <- regions
  if(!missing(sites)) {
    res$Represented <- as.character(res[,col]) %in% as.character(unique(sites[,col]))
    N_studies <- tapply(sites$SS, 
                        sites[,col], 
                        function(ss) length(unique(ss)))
    stopifnot(all(names(N_studies) %in% res[,col]))
    res$N_studies <- N_studies[as.character(res[,col])]
    # Empty names are not matched
    if('' %in% sites[,col]) {
      res$N_studies[''==res[,col]] <- N_studies[''==names(N_studies)]
    }
    res$N_studies[is.na(res$N_studies)] <- 0
    res$Percent_studies <- 100*res$N_studies/sum(res$N_studies)
    stopifnot(all.equal(sum(res$Percent_studies), 100))

    N_sites <- table(sites[,col])
    stopifnot(all(names(N_sites) %in% res[,col]))
    res$N_sites <- as.numeric(N_sites[as.character(res[,col])])
    res$N_sites[is.na(res$N_sites)] <- 0
    # Empty names are not matched
    if('' %in% sites[,col]) {
      res$N_sites[''==res[,col]] <- N_sites[''==names(N_sites)]
    }
    res$Percent_sites <- 100*res$N_sites/sum(res$N_sites, na.rm=TRUE)
    stopifnot(all.equal(sum(res$Percent_studies), 100))
  }

  if(!missing(diversity)) {
    N_samples <- table(diversity[,col])
    stopifnot(all(names(N_samples) %in% res[,col]))
    res$N_samples <- as.numeric(N_samples[as.character(res[,col])])
    # Empty names are not matched
    if('' %in% diversity[,col]) {
      res$N_samples[''==res[,col]] <- N_samples[''==names(N_samples)]
    }
    res$N_samples[is.na(res$N_samples)] <- 0
    res$Percent_samples <- 100*res$N_samples/sum(res$N_samples, na.rm=TRUE)
    stopifnot(all.equal(sum(res$Percent_samples, na.rm=TRUE), 100))
  }
  return (res)
}

RealmsTable <- function(sites, diversity) {
  # Takes site-level and diversity-level database extracts and returns a 
  # data.frame of realm-level statistics.
  data(realms)
  return (.RegionsTable(realms, 'Realm', sites, diversity))
}

BiomesTable <- function(sites, diversity) {
  # Takes site-level and diversity-level database extracts and returns a 
  # data.frame of biome-level statistics.
  data(biomes)
  return (.RegionsTable(biomes, 'Biome', sites, diversity))
}

EcoregionsTable <- function(sites, diversity) {
  # Takes site-level and diversity-level database extracts and returns a 
  # data.frame of eco-region-level statistics.
  data(ecoregions)
  return (.RegionsTable(ecoregions, 'Ecoregion', sites, diversity))
}

CountriesTable <- function(sites, diversity) {
  # Takes site-level and diversity-level database extracts and returns a 
  # data.frame of country-level statistics.
  data(countries)
  return (.RegionsTable(countries, 'Country', sites, diversity))
}

HotspotsTable <- function(sites, diversity) {
  # Takes a site-level and diversity database extracts and returns a data.frame
  # of hotspot-level statistics.
  data(hotspots)
  percent.non.hotspot <- 100-sum(hotspots$Percent_terrestrial_area)
  h <- rbind(data.frame(Hotspot='',
                        Area_km2=NA,
                        Percent_terrestrial_area=percent.non.hotspot,
                        Predominant_realm=''),
             hotspots)
  hotspots$Predominant_realm <- EnsureFactorLevels(hotspots$Predominant_realm,
      Realms(), allow.empty=TRUE)
  return (.RegionsTable(h, 'Hotspot', sites, diversity))
}

WildernessAreasTable <- function(sites, diversity) {
  # Takes a site-level and diversity database extracts and returns a data.frame
  # of wilderness area-level statistics.
  data(wilderness.areas)
  percent.non.wilderness.areas <- 100 - sum(wilderness.areas$Percent_terrestrial_area)
  wa <- rbind(data.frame(Wilderness_area='',
                         Area_km2=NA,
                         Perimeter=NA,
                         Percent_terrestrial_area=percent.non.wilderness.areas,
                         Predominant_realm=''),
              wilderness.areas)
  return (.RegionsTable(wa, 'Wilderness_area', sites, diversity))
}

PlotBiomesCoverage <- function(sites, diversity, mar=c(3.5,4.3,0.8,0.5), 
                               mgp=c(2.2,0.6,0), las=1, ...) {
  op <- par(mar=mar, mgp=mgp, mfrow=c(3,2), las=las)
  on.exit(par(op))

  Plot <- function(rows, xcol, ycol, xlab, ylab, col, labels) {
      # Helper
      x <- log10(rows[,xcol])
      y <- log10(rows[,ycol])
      # Increase axis limits to accomodate the filled box under letters
      xlim <- c(min(x) - diff(range(x)) * 0.03, max(x) + diff(range(x)) * 0.02)
      ylim <- c(min(y) - diff(range(y)) * 0.03, max(y) + diff(range(y)) * 0.02)
      plot(x, y, xlab=xlab, ylab=ylab, type='n', xlim=xlim, ylim=ylim, ...)
      OneToOneLine()
      points(x, y, col=as.character(rows$Weak_colour), pch=15, cex=3, ...)
      text(x, y, rows$ID, ...)
      sub.heading <<- SubHeading(sub.heading)

  }
  OneToOneLine <- function() abline(a=0,  b=1, col='grey', ...)

  sub.heading <- 1
  biomes <- BiomesTable(sites, diversity)

  Plot(biomes[biomes$Percent_studies>0,],
       'Percent_terrestrial_NPP', 'Percent_studies',
       xlab='', ylab=~log[10]('%'~Studies))
  Plot(biomes[biomes$Percent_studies>0,],
       'Percent_terrestrial_area', 'Percent_studies',
       xlab='', ylab='')

  Plot(biomes[biomes$Percent_sites>0,],
       'Percent_terrestrial_NPP', 'Percent_sites',
       xlab='', ylab=~log[10]('%'~Sites))
  Plot(biomes[biomes$Percent_sites>0,],
       'Percent_terrestrial_area', 'Percent_sites',
       xlab='', ylab='')

  Plot(biomes[biomes$Percent_samples>0,],
       'Percent_terrestrial_NPP', 'Percent_samples',
       xlab=~log[10]('%'~NPP), ylab=~log[10]('%'~Samples))
  Plot(biomes[biomes$Percent_samples>0,],
       'Percent_terrestrial_area', 'Percent_samples',
       xlab=~log[10]('%'~area), ylab='')
  par(op)
}

.PlotCountryPoints <- function(values, x, y, mar, mgp, cex, labels.cex, 
                               show.labels, ...) {
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))

  require(wordcloud)

  values <- values[values[,x]>0 & values[,y]>0,]

  # Colours
  col <- rep('white', nrow(values)) # Countries not in tropics are white
  col.pool <- rev(heat.colors(10))
  tropical <- values$Fraction_area_in_tropics>0
  col[tropical] <- col.pool[ceiling(10*values[tropical,'Fraction_area_in_tropics'])]

  x <- as.vector(log10(values[,x]))
  y <- as.vector(log10(values[,y]))

  # Empty plot
  plot(x, y, type='n', ...)
  axis(side=3, labels=FALSE, ...)
  axis(side=4, labels=FALSE, ...)

  # Label positions
  nc <- wordlayout(c(x,x), c(y,y), c(rep('X', nrow(values)), values$ISO_2_digit),
                   xlim=range(x), ylim=range(y), cex=labels.cex)
  nc <- nc[(1+(nrow(nc)/2)):nrow(nc),]

  # Lines
  if(show.labels) segments(x, y, nc[,1]+nc[,3]/2, nc[,2]+nc[,4]/2, col='lightgrey', ...)

  # Points
  points(x, y, pch=ifelse(values$Megadiverse, 24, 21), bg=col, col=1, cex=cex, ...)

  # Labels
  if(show.labels) text(nc[,1] + .5*nc[,3],nc[,2]+.5*nc[,4],values$ISO_2_digit, cex=labels.cex, ...)
}

PlotNStudiesAgainstCountryArea <- function(sites, 
                                           mar=c(4,4.3,0.5,0.5), 
                                           mgp=c(2.5, 0.6, 0), 
                                           cex=1, labels.cex=0.7, 
                                           xlab=~log[10](country~area~(km^2)), 
                                           ylab=~log[10]('%'~Studies), 
                                           show.labels=TRUE, 
                                           ...) {
  .PlotCountryPoints(CountriesTable(sites=sites), 'Area', 'Percent_studies', 
                     mar, mgp, cex, labels.cex, show.labels, 
                     xlab=xlab, ylab=ylab, ...)
}

PlotNSitesAgainstCountryArea <- function(sites, 
                                         mar=c(4,4.3,0.5,0.5), 
                                         mgp=c(2.5, 0.6, 0), 
                                         cex=1, labels.cex=0.7, 
                                         xlab=~log[10](country~area~(km^2)), 
                                         ylab=~log[10]('%'~Sites), 
                                         show.labels=TRUE,
                                         ...) {
  .PlotCountryPoints(CountriesTable(sites=sites), 'Area', 'Percent_sites', 
                     mar, mgp, cex, labels.cex, show.labels, 
                     xlab=xlab, ylab=ylab, ...)
}

PlotNSamplesAgainstCountryArea <- function(diversity, 
                                           mar=c(4,4.3,0.5,0.5), 
                                           mgp=c(2.5, 0.6, 0), 
                                           cex=1, labels.cex=0.7, 
                                           xlab=~log[10](country~area~(km^2)), 
                                           ylab=~log[10]('%'~Samples), 
                                           show.labels=TRUE, 
                                           ...) {
  .PlotCountryPoints(CountriesTable(diversity=diversity), 'Area', 
                     'Percent_samples', mar, mgp, cex, labels.cex, show.labels, 
                     xlab=xlab, ylab=ylab, ...)
}

.BinBars <- function(breaks, border=NA, col='lightgrey', ...) {
  # Plots bars showing bins
  usr <- par('usr')
  rect(usr[1], 
       head(breaks[seq(1, length(breaks), 2)], -1), 
       usr[2], 
       breaks[seq(2, length(breaks), 2)], 
       border=border, col=col, ...)
}

PlotLatitudes <- function(sites, diversity, sanity=NULL, show.bin.mids=FALSE, 
                           show.bin.bars=FALSE, ...) {
  mgp <- c(2.5,0.75,0)
  op <- par(mfrow=c(1,2), 
            mar=c(4,0.4,0.5,2),     # bottom, left, top, right
            mgp=mgp)
  on.exit(par(op))

  # Compute N sites, studies and samples in 5-degree bands
  data(terrestrial.bands)
  stopifnot(all(5==diff(terrestrial.bands$Lat_lower)))
  breaks <- seq(-90, 90, by=5)
  y.axis.ticks <- seq(-90, 90, by=10)

  # Compute N sites, studies and samples per band
  sites_ <- hist(sites$Latitude, breaks=breaks, plot=FALSE)
  sites_$count <- 100*sites_$count/sum(sites_$count)
  sites_$count[0==sites_$count] <- NA

  studies <- hist(tapply(sites$Latitude, sites$SS, median), breaks=sites_$breaks,
                  plot=FALSE)
  studies$count <- 100*studies$count/sum(studies$count)
  studies$count[0==studies$count] <- NA

  samples <- hist(diversity$Latitude, breaks=sites_$breaks, plot=FALSE)
  samples$count <- 100*samples$count/sum(samples$count)
  samples$count[0==samples$count] <- NA

  # Extend y axis to show all sites, studies and samples
  ylim <- range(sites_$mids[!is.na(sites_$count) | !is.na(studies$count) | !is.na(samples$count)])

  plot(sites_$count, sites_$mids, type='n', 
       xlim=rev(range(c(sites_$count, studies$count, samples$count), na.rm=TRUE)), 
       ylim=ylim, xlab='% studies/sites/samples', ylab='', axes=FALSE, ...)

  if(show.bin.bars) {
    .BinBars(sites_$breaks)
  }

  # axis lines
  usr <- par('usr')
  lines(c(usr[2], usr[2]), c(usr[3], usr[4]), xpd=TRUE)
  lines(c(usr[1], usr[2]), c(usr[3], usr[3]), xpd=TRUE)
  axis(side=1, ...)
  axis(side=4, at=y.axis.ticks, labels=FALSE, ...)


  if(show.bin.mids) {
    # Bin mid-points
    abline(h=studies$mids, col='lightgrey')
  }

  # Tropics
  abline(h=c(-23.5,23.5), col='darkgrey', lty=3)

  if(!is.null(sanity)) {
    abline(h=sanity)
    abline(v=sites_$count[sanity==sites_$mids])
    abline(v=studies$count[sanity==sites_$mids])
    abline(v=samples$count[sanity==sites_$mids])
  }

  points(studies$count, studies$mids, pch=19, ...)
  points(sites_$count, sites_$mids, pch=4, ...)
  points(samples$count, samples$mids, pch=3, ...)

  if(!is.null(sanity)) abline(v=studies$counts[sanity==studies$mid])

  # Terrestrial area by latitudinal band
  par(mar=c(4,1,0.5,0.5), mgp=mgp)

  with(terrestrial.bands, {
    xlim <- range(c(Percent_terrestrial_area, Percent_global_terrestrial_NPP))
    plot(NA, NA, xlab='% total terrestrial area/NPP', ylab='', axes=FALSE, 
         lty=3, ylim=ylim, xlim=xlim, ...)

    if(show.bin.bars) {
      .BinBars(sites_$breaks)
    }

    # axis lines
    usr <- par('usr')
    lines(c(usr[1], usr[1]), c(usr[3], usr[4]), xpd=TRUE)
    lines(c(usr[1], usr[2]), c(usr[3], usr[3]), xpd=TRUE)
    axis(side=1, ...)
    axis(side=2, at=y.axis.ticks, ...)

    if(show.bin.mids) {
      # Bin mid-points
      abline(h=studies$mids, col='lightgrey')
    }

    lines(head(rep(Percent_terrestrial_area, each=2), -1), 
          tail(rep(Lat_lower, each=2),-1))

    lines(head(rep(Percent_global_terrestrial_NPP, each=2), -1), 
          tail(rep(Lat_lower, each=2),-1), lty=2)

    if(!is.null(sanity)) {
      abline(h=sanity)
      abline(v=Percent_terrestrial_area[sanity==Lat_centre])
      abline(v=Percent_global_terrestrial_NPP[sanity==Lat_centre])
    }
  })

  # Tropics
  abline(h=c(-23.5,23.5), col='darkgrey', lty=3)

  return (invisible(sites_))
}

PlotSamplingDurations <- function(sites) {
  op <- par(mar=c(3,3,0.5,0.5), mgp=c(2,0.3,0))
  on.exit(par(op))

  duration <- log10(as.numeric(sites$Sample_end_latest - sites$Sample_start_earliest))
  n.infinite <- sum(is.infinite(duration))
  duration <- duration[!is.infinite(duration)]
  hist(duration,xlab=~log[10](days), main='')
  mean_ <- mean(duration)
  median_ <- median(duration)
  abline(v=c(mean_, median_), col='darkgrey', lty=1:2)
  return (invisible(c(n.infinite=n.infinite, mean=10^mean_, median=10^median_)))
}

PlotStudyCommonTaxa <- function(taxa, mar=c(2,11.5,0.5,8), mgp=c(1.4,0.3,0),
                                cex.names=0.7, cex.axis=0.9, border, xaxs='i', 
                                yaxs='i', show.references=TRUE, ...) {
  # xaxs and yaxs default to 'i' (internal) to prevent the gap of 4% of the 
  # data's range from being added between the data and the axis. See help page 
  # for par.
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))

  # Subset taxa
  ts <- unique(taxa[,c('Source_ID','SS','Study_common_taxon',
                       'Rank_of_study_common_taxon')])
  stopifnot(nrow(ts)==length(unique(ts$SS)))

  # Show ranks in descending order
  ts$Rank_of_study_common_taxon <- factor(ts$Rank_of_study_common_taxon, 
                                          levels=c('',rev(TaxonomicRanks())))
  ts$Rank_of_study_common_taxon <- droplevels(ts$Rank_of_study_common_taxon)

  # N studies for that group
  res <- data.frame(Rank=tapply(ts$Rank_of_study_common_taxon, ts$Study_common_taxon, unique),
                    N=as.vector(table(ts$Study_common_taxon)))

  stopifnot(sum(res[,'N'])==length(unique(ts$SS)))

  # Order by decreasing rank, alphabetically where tied
  res <- res[order(res[,'Rank'], -res[,'N'], rownames(res), decreasing=TRUE),]

  # Colours
  colours <- TaxonomicGroupSoftColours()
  col <- colours[as.character(taxa$Kingdom[match(rownames(res), taxa$Study_common_taxon)])]
  col[is.na(col)] <- colours[as.character(taxa$Phylum[match(rownames(res)[is.na(col)], taxa$Study_common_taxon)])]
  col[is.na(col)] <- colours['Other']
  col['Animalia'==rownames(res)] <- colours['Other']

  # Empty indicates more than one kingdom
  col[''==rownames(res)] <- 'black'
  rownames(res)[''==rownames(res)] <- 'Multiple kingdoms'

  if(missing(border))     border <- col

  y <- HorizontalBarplot(res[,'N'], show.labels=FALSE, cex.names=cex.names, 
                         cex.axis=cex.axis, xlab=~Number~of~Studies,
                         yaxt='n', col=col, border=border, 
                         xaxs=xaxs, yaxs=yaxs, ...)

  # Bands indicating common rank
  bands <- tapply(y, res[,'Rank'], range)
  bands <- do.call('rbind', bands)

  # Plot outside of the normal plotting area
  par(xpd=NA)

  # X axis line across the width of the plot
  segments(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[3])

  rect(grconvertX(0, 'ndc', 'user'), 
       bands[,1]-diff(y)[1]/2, 
       -0.15,   # Nasty hack to avoid over-printing bars - s
       bands[,2]+diff(y)[1]/2, 
       col=c('white','lightgrey'), border=NA, ...)

  # Taxon names
  labels <- gsub(' ', '~', rownames(res))
  rank <- levels(ts$Rank_of_study_common_taxon)[res[,'Rank']]
  italics <- rank %in% c('Genus','Species','Infraspecies')
  labels[italics] <- paste('italic(', labels[italics], ')', sep='')
  labels[''==labels] <- '""'
  labels[''!=labels] <- paste('~', labels[''!=labels], sep='')
  # Should be possible to do this using axis but I couldn't get axis to 
  # work with expressions
  labels <- sapply(sapply(sapply(labels, as.expression), eval), as.formula)
  mapply(FUN=mtext, labels, at=y, side=2, adj=1, las=1, line=0, 
         cex=cex.names, ...)

  # Group names
  labels <- rev(levels(ts$Rank_of_study_common_taxon)[unique(res[,'Rank'])])
  n <- paste('(', tapply(res[,'N'], res[,'Rank'], sum), ')', sep='')  
  labels[bands[,1]!=bands[,2]] <- paste(labels, n)[bands[,1]!=bands[,2]]

  text(grconvertX(0, 'ndc', 'user'), apply(bands, 1, mean), labels, 
       adj=c(-0.05,NA), cex=cex.axis, ...)
  par(xpd=FALSE)

  # Assign reference numbers
  rownames(res)['Multiple kingdoms'==rownames(res)] <- ''
  references <- lapply(rownames(res), function(t) {
      data.frame(Taxon=t, 
                 Source_ID=as.character(unique(taxa$Source_ID[t==taxa$Study_common_taxon])))
  })

  references <- do.call('rbind', references)
  references <- references[nrow(references):1,]

  # Assign numerical ID based on the author name and year in the source ID
  reference.id <- paste0(substring(references$Source_ID, 11),      # First author
                         substring(references$Source_ID, 5, 8),    # Year
                         references$Source_ID)                     # Source ID to 
                                                                   # avoid collisions
  references <- cbind(references, ID=reference.id)

  # Order by author, year within common taxon
  references <- references[order(-as.integer(references$Taxon), references$ID),]
  references$ID <- as.integer(factor(references$ID, levels=unique(references$ID)))

  if(show.references) {
    # Lines the ends of the bars to the reference numbers
    centres <- cbind(res[,'N']+0.2, y[,1])
    centres <- centres[seq(1,nrow(res),by=2),,drop=FALSE]
    segments(centres[,1], centres[,2], par('usr')[2], centres[,2], lty=1, col='grey')

    Group <- function(values) {
      values <- sort(values)
      last <- head(values,1)
      groups <- 1
      for(value in tail(values,-1)) {
        group <- tail(groups,1)
        if(value-last>1) group <- group+1
        groups <- c(groups, group)
        last <- value
      }

      tapply(values, groups, function(values) {
        ifelse(1==length(unique(values)), values, 
                                  paste(head(values,1),'-',tail(values,1),sep=''))
      })
    }

    labels <- tapply(references$ID, references$Taxon, Group)
    labels <- sapply(labels, paste, collapse=',')
    axis(labels, side=4, at=y, las=1, tick=FALSE, cex.axis=cex.names, ...)
  }

  studies <- tapply(ts$SS, ts$Study_common_taxon, unique)
  sources <- tapply(ts$Source_ID, ts$Study_common_taxon, unique)
  invisible(list(y=y, references=references, studies=studies, sources=sources))
}


PlotMatrix <- function(x, lim=range(x), mar=c(0.5,9,2.5,0.5), mgp=c(3, 0.7, 0), 
                       cex.axis=par('cex.axis'), colour.ramp=heat.colors(10), 
                       labels=NULL, ...) {
  # Plot a matrix as a coloured image
  # Inspired by http://www.phaget4.org/R/myImagePlot.R
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))

  ylab <- rownames(x)
  xlab <- colnames(x)

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

  colour.levels <- seq(lim[1], lim[2], length=length(colour.ramp))

  # Reverse Y axis
  reverse <- nrow(x):1
  ylab <- ylab[reverse]
  x <- x[reverse,]

  # Data Map
  # Expand zlim to match colour scale
  bin.width <- ((lim[2]-lim[1]) / (length(colour.ramp)-1))/2
  image(1:length(xlab), 1:length(ylab), t(x), col=colour.ramp, xlab="",
        ylab="", axes=FALSE, zlim=c(lim[1]-bin.width, lim[2]+bin.width), ...)
  box()
  axis(side=3, at=1:length(xlab), labels=xlab, las=1, cex=cex.axis, pos=NA, 
       tick=FALSE, ...)
  axis(side=2, at=1:length(ylab), labels=ylab, las=1, cex=cex.axis, pos=NA, 
       tick=FALSE, ...)

  if(!is.null(labels)) {
    text(x=col(x), y=row(x), labels[reverse,], ...)
  }

  # Colour Scale
  par(mar=c(mar[1], 4, mar[3:4]))
  image(1, colour.levels,
        matrix(data=colour.levels, ncol=length(colour.levels), nrow=1),
        col=colour.ramp, xlab="", ylab="", xaxt="n", yaxt='n')
  mtext('Difference', side=3, line=par('mgp')[2])
  axis(side=2, at=colour.levels, las=1, cex.axis=cex.axis, ...)
  box()

  layout(1)
}

PlotLULC <- function(sites, lulc.2005, lim=c(-20,20),
                     colour.ramp=c('#2166ac','#4393c3','#92c5de','#d1e5f0',
                                  '#ffffff','#fddbc7','#f4a582','#d6604d',
                                  '#b2182b'), ...) {
  # An image comparing the land-use and land-cover distributions of sites with 
  # estimate for 2005.

  sites <- table(sites$Predominant_habitat, sites$Use_intensity)
  stopifnot(all(2:5 == grep('*secondary*', rownames(sites), ignore.case=TRUE)))
  stopifnot('Plantation forest'==rownames(sites[6]))
  sites <- rbind(sites[1,,drop=FALSE], 
                 'Secondary vegetation'=colSums(sites[2:5,]), 
                 sites[6:nrow(sites),])

  intensities <- c('Minimal use', 'Light use', 'Intense use')
  habitats <- c('Primary vegetation','Secondary vegetation','Plantation forest',
                'Cropland','Pasture','Urban')
  sites <- 100*sites/sum(sites)

  # Tim placed all plantation forest in minimal - spread plantation forest values 
  # evenly
  stopifnot(all(is.na(lulc.2005['plantation forest',c('light','intense')])))
  lulc.2005['plantation forest',] <- lulc.2005['plantation forest','minimal']/3

  # Adjust names
  stopifnot(all(rownames(lulc.2005)[1:2]==c('primary','secondary')))
  rownames(lulc.2005)[1:2] <- paste(rownames(lulc.2005)[1:2], 'vegetation')
  colnames(lulc.2005) <- paste(Capitalize(colnames(lulc.2005)), 'use')
  rownames(lulc.2005) <- Capitalize(rownames(lulc.2005))

  # Percentages
  lulc.2005 <- 100*lulc.2005/sum(lulc.2005)

  # Should be percentages
  stopifnot(100==sum(sites))
  stopifnot(100==sum(lulc.2005))

  # Should have same col and row names
  sites <- sites[habitats,intensities]
  lulc.2005 <- lulc.2005[habitats,intensities]

  difference <- as.matrix(sites-lulc.2005)
  difference.str <- FormatPercent(difference)
  difference.str[0<difference] <- paste('+', difference.str[0<difference], sep='')
  difference.str <- paste0('(', difference.str, ')')

  # Don't show difference where we don't have an estimate
  difference.str[0==lulc.2005] <- ''

  labels <- paste0(FormatPercent(as.matrix(sites)), '\n', difference.str)
  dim(labels) <- dim(difference)
  PlotMatrix(difference, lim=lim, colour.ramp=colour.ramp, labels=labels, ...)
  invisible (difference)
}

PlotJournals <- function(bib, mar=c(2,11.5,0.5,0.5), mgp=c(1,0.3,0), 
                         cex.names=0.7, cex.axis=0.7, yaxs='i',
                         xlab='Percentage', ...) {
  # yaxs defaults to 'i' (internal) to prevent the gap of 4% of the 
  # data's range from being added between the data and the axis. See help page 
  # for par.
  op <- par(mar=mar, mgp=mgp)
  on.exit(par(op))

  journals <- cbind(Sources=table(bib$Journal_title), 
                    Studies=tapply(bib$N_studies, bib$Journal_title, sum), 
                    Sites=tapply(bib$N_sites, bib$Journal_title, sum), 
                    Samples=tapply(bib$N_samples, bib$Journal_title, sum))
  journals <- 100*sweep(journals, 2, colSums(journals), '/')

  # Order by n sources
  journals <- journals[order(journals[,'Sources']),]
  journals <- journals[journals[,'Sources']>1 & ''!=rownames(journals),]
  y <- HorizontalBarplot(journals[,'Sources'], show.labels=FALSE, 
                         cex.names=cex.names, cex.axis=cex.axis, xlab=xlab,
                         xlim=range(journals), yaxs=yaxs, xaxs='r', ...)
  points(journals[,'Studies'], y, pch=19)
  points(journals[,'Sites'], y, pch=4)
  points(journals[,'Samples'], y, pch=3)
}

Histogram <- function(v, mar=c(4,4,2,0.5), mgp=c(3, 0.7, 0), 
                      xaxs='i', yaxs='i', N.cex=2, las=1, ...) {
  # A slightly nicer histogram than the usual
  op <- par(mar=mar, mgp=mgp, xaxs=xaxs, yaxs=yaxs, las=las)
  on.exit(par(op))
  res <- hist(v, ...)

  m <- median(v, na.rm=TRUE)
  abline(v=m)
  axis(side=3, at=m, labels=sprintf('median=%.2f', m), ...)

  text(par('usr')[2], par('usr')[4], 
       paste('N=', FormatThousands(sum(!is.na(v) & !is.infinite(v))), sep=''), 
       cex=N.cex, adj=c(1.2,1.2), ...)

  invisible(res)
}

PlotNSitesByCollationDate <- function(bib, mar=c(4,4,1.5,0.5), mgp=c(3, 0.7, 0),
                                      las=1, ...) {
    # Number of sites by date entered into the database.
    op <- par(mar=mar, mgp=mgp, las=las)
    on.exit(op)
    events <- data.frame(
        Label=c('Projects start','Projects end','Projects start','Projects end','Interns start','Interns end','Projects start'),
        Date=as.Date(c('2012-02-01','2012-09-01','2013-02-01','2013-09-01','2013-10-10','2014-04-01','2014-06-01')),
        Line=c(0, 0.5, 0, 0, 0.5, 0, 0.5)
    )
    with(bib[order(bib$Compiled_on),],
    {
        plot(Compiled_on, cumsum(N_sites), xlab='Date',
             xlim=range(c(events$Data, Compiled_on)), ylab='N sites', pch=19, xaxt='n')
    })

    axis.Date(
        1, at=seq(
            ISOdate(substr(min(c(events$Date, bib$Compiled_on)), 1, 4), 1, 1),
            ISOdate(1+as.numeric(substr(max(c(events$Date, bib$Compiled_on)), 1, 4)), 1, 1),
            by="6 months"),
        format='%b %Y',
        ...
    )

    abline(v=events$Date)
    mtext(events$Label, side=3, at=events$Date, line=events$Line, cex=0.8)
}

PlotBiomeByDate <- function(sites, diversity, mar=c(3, 3, 0.5, 0.5),
                            mgp=c(2, 0.5, 0), las=1, col='#00000010',
                            x.lim.lower.padding=365*2,
                            coverage.label.offset=500,
                            ...) {
    # x.lim.lower.padding: xlim[1] is computed as min(sites$Sample_start_earliest) - x.lim.lower.padding,
    # to allow space for the text showing n samples and n sites
    # coverage.label.offset: text showing n samples and n sites per biome is plotted at xlim[1] + coverage.label.offset
    op <- par(mar=mar, mgp=mgp, las=las)
    on.exit(par(op))

    data(biomes)

    n.samples <- table(diversity$SSS)
    sites$N_samples <- n.samples[match(sites$SSS, names(n.samples))]
    sites$N_samples[is.na(sites$N_samples)] <- 0
    sites$log10_N_samples <- log10(1+sites$N_samples)
    sites$Sample_date_mid <- sites$Sample_start_earliest+(sites$Sample_end_latest-sites$Sample_start_earliest)/2

    # Show biome IDs on y axis
    ids <- as.character(sort(biomes$ID, decreasing=TRUE))

    # Don't show 'Inland Water' or 'Rock & Ice'
    stopifnot(all(0==table(sites$Biome)[c('Inland Water','Rock & Ice')]))
    ids <- tail(ids, -2)
    sites$Biome <- as.character(sites$Biome)
    sites$Biome <- factor(sites$Biome, levels=setdiff(Biomes(), c('Inland Water','Rock & Ice')))
    diversity$Biome <- as.character(diversity$Biome)
    diversity$Biome <- factor(diversity$Biome, levels=setdiff(Biomes(), c('Inland Water','Rock & Ice')))

    BiomeNameToY <- function(biome) {
        # Match site-level biome name to row in biomes data.frame
        biome.index <- match(biome, biomes$Biome)
        return (match(biomes$ID[biome.index], ids))
    }

    SanitizeBiomeName <- function(n) {
        n <- gsub(', ', ',\n', n)
        n <- gsub('Subtropical ', 'Subtropical\n', n)
    }

    # Compute limits of x axis
    from.date <- min(sites$Sample_start_earliest)
    # Margin on left for N sites and N samples
    from.date <- from.date - x.lim.lower.padding
    to.date <- max(sites$Sample_end_latest)

    plot(sites$Sample_date_mid, as.numeric(sites$Biome), type='n',
         xlim=c(from.date, to.date), xlab='Year', xaxt='n',
         ylim=c(1, nlevels(sites$Biome)), ylab='', yaxt='n',
         frame.plot=FALSE, ...)

    x.ticks <- seq(as.Date(ISOdate(substr(from.date, 1, 4), 1, 1)),
                   as.Date(ISOdate(1 + as.numeric(substr(to.date, 1, 4)), 1, 1)),
                   by="year")
    axis(1, at=as.numeric(as.Date(x.ticks)), labels=FALSE, lwd=0, lwd.ticks=1, ...)

    x.labels <- seq(as.Date(ISOdate(substr(from.date, 1, 4), 7, 1)),
                    as.Date(ISOdate(substr(to.date, 1, 4), 7, 1)),
                    by="year")
    axis(side=1, at=as.numeric(as.Date(x.labels)), labels=format(x.labels, '%Y'),
         cex.axis=0.7, tick=FALSE, line=-0.5, ...)

    # Previous solution for reference
    # Rotate x axis labels
    # http://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-I-create-rotated-axis-labels_003f
    # text(x=as.numeric(as.Date(x)), y=par("usr")[3] - 0.3,
    #      labels=format(x, '%Y'), xpd=NA, srt=35, cex=0.7, ...)

    # x axis title at top
    mtext('Biome', side=2, at=par('usr')[4], line=5)

    axis(side=2, at=1:nlevels(sites$Biome), labels=ids, las=1, line=10.5,
         tick=FALSE, cex.axis=1, ...)
    axis(side=2, at=1:nlevels(sites$Biome), las=1, line=-0.25, tick=FALSE,
         cex.axis=0.75,
         labels=SanitizeBiomeName(as.character(biomes$Biome)[match(ids, biomes$ID)]),
         ...)

    y <- (1:nlevels(sites$Biome))
    rect(par('usr')[1], y-0.5, par('usr')[2], y+0.5, border=NA,
          col=as.character(biomes$Weak_colour)[match(ids, biomes$ID)], ...)

    y <- BiomeNameToY(sites$Biome)
    jitter <- rnorm(n=nlevels(sites$SS), sd=0.1)
    y <- y + jitter[as.numeric(sites$SS)]

    segments(sites$Sample_start_earliest, y, sites$Sample_end_latest, y,
             col=col, ...)

    symbols(sites$Sample_date_mid, y, circles=60*sqrt(sites$log10_N_samples),
            bg=NA, fg=col, add=TRUE, inches=FALSE)

    text(from.date + coverage.label.offset,
         BiomeNameToY(levels(sites$Biome)),
         paste('', FormatThousands(table(sites$Biome)), 'sites\n',
               FormatThousands(table(diversity$Biome)), 'samples'),
         adj=c(1,0.5), cex=0.7, ...)
}

PlotLatitudeByDate <- function(sites, diversity, mar=c(3,3,0.5,0.5),
                               mgp=c(2,0.5,0), las=1, ...) {
    op <- par(mar=mar, mgp=mgp, las=las)
    on.exit(par(op))

    n.samples <- table(diversity$SSS)
    sites$N_samples <- n.samples[match(sites$SSS, names(n.samples))]
    sites$N_samples[is.na(sites$N_samples)] <- 0
    sites$log10_N_samples <- log10(1+sites$N_samples)
    sites$Sample_date_mid <- sites$Sample_start_earliest+(sites$Sample_end_latest-sites$Sample_start_earliest)/2

    from.date <- min(sites$Sample_start_earliest)
    to.date <- max(sites$Sample_end_latest)

    plot(sites$Sample_date_mid, abs(sites$Latitude), type='n',
         xlim=c(from.date, to.date), xlab='', xaxt='n', ylab=~abs(latitude), ...)
    axis(1, at=seq(from.date, 365+to.date, by="year"), labels=FALSE, ...)

    # Rotate x axis labels
    # http://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-I-create-rotated-axis-labels_003f
    x <- seq(ISOdate(substr(from.date, 1, 4), 6, 1),
             ISOdate(substr(to.date, 1, 4), 6, 1), by="year")
    text(x=as.numeric(as.Date(x)), y=par("usr")[3]-4,
                      labels=format(x, '%Y'), xpd=NA, srt=35, cex=0.7, ...)

    data(realms)

    col <- AddAlpha(as.character(realms$Colour), 0x10)
    col <- col[match(sites$Realm, realms$Realm)]
    segments(sites$Sample_start_earliest, abs(sites$Latitude),
             sites$Sample_end_latest, abs(sites$Latitude), col=col, ...)

    if(FALSE) {
      # The old way, plotting points
      points(sites$Sample_date_mid, abs(sites$Latitude),
             cex=sites$log10_N_samples*0.8, col=col, bg=col, pch=21, ...)
    } else {
      symbols(sites$Sample_date_mid, abs(sites$Latitude),
              circles=60*sqrt(sites$log10_N_samples), bg=col, fg=NA, add=TRUE,
              inches=FALSE)
    }

    stopifnot(!any('Antarctic'==sites$Realm))
    with(droplevels(realms['Antarctic'!=realms$Realm,]),
         legend("topleft", legend=levels(Realm), col=as.character(Colour),
         pch=19, cex=0.7, pt.cex=1, y.intersp=0.7))
}

PlotSpatioTemporal <- function(sites, left.mar=4, bottom.mar=4.05,
                               absolute.lat=FALSE, las=1, sanity.lat=NULL,
                               sanity.date=NULL, ...) {
    #                  left   right   bottom  top
    screen <- matrix(c(0,     3/4,    0,      3/4,
                       0,     3/4,    3/4,    1,
                       3/4,   1,      0,      3/4),
                     ncol=4, byrow=TRUE)

    split.screen(screen)
    screen(1)

    op <- par(mar=c(bottom.mar, left.mar, 0, 0), las=las)
    on.exit(par(op))

    sites$Sample_date_mid <- sites$Sample_start_earliest+(sites$Sample_end_latest-sites$Sample_start_earliest)/2

    if(absolute.lat) {
        lat <- abs(sites$Latitude)
        lat.lab <- ~abs(latitude)
    } else {
        lat <- sites$Latitude
        lat.lab <- 'latitude'
    }

    # Compute histograms to get axes ranges
    date.h <- hist(sites$Sample_date_mid, breaks='years', plot=FALSE)
    lat.h <- hist(lat, plot=FALSE)

    date.range <- range(date.h$breaks)
    lat.range <- range(lat.h$breaks)
    plot(sites$Sample_date_mid, lat, pch=19, xlim=date.range,
         ylim=lat.range, xlab='Sampling date mid-point', ylab=lat.lab)
    if(!is.null(sanity.lat)) {
        abline(h=sanity.lat)
    }
    if(!is.null(sanity.date)) {
        abline(v=sanity.date)
    }

    screen(2)
    par(mar=c(0, left.mar, 0, 0))
    plot(range(sites$Sample_date_mid), range(date.h$counts), type='n', xaxt='n',
         yaxt='s', main='', frame=FALSE, ylab='Count', xlim=date.range)
    rect(tail(date.h$breaks, -1), 0, head(date.h$breaks, -1), date.h$counts)

    screen(3)
    par(mar=c(bottom.mar, 0, 0, 0))
    plot(range(lat.h$counts), lat.range, type='n', xaxt='s', yaxt='n', main='',
         frame=FALSE, xlab='Count')
    rect(0, tail(lat.h$breaks, -1), lat.h$counts, head(lat.h$breaks, -1))

    close.screen(all=TRUE)
}

PlotHotspotsCoverage <- function(sites, diversity, mar=c(4, 20, 0.5, 0.5),
         mgp=c(2, 0.6, 0),         # axis title, axis labels and axis line
         ...) {
    # A horizontal barplot of hotspots coverage
    op <- par(mar=mar, mgp=mgp)
    on.exit(par(op))

    hotspots <- HotspotsTable(sites, diversity)
    hotspots <- hotspots[''!=hotspots$Hotspot,]    # Drop 'no hotspot'
    hotspots$Area_corrected_sites = hotspots$N_sites / hotspots$Area_km2
    hotspots$Area_corrected_samples = hotspots$N_samples / hotspots$Area_km2
    hotspots <- hotspots[nrow(hotspots):1,]

    hotspots$Percent_hotspot_area <- 100 * hotspots$Area_km2 / sum(hotspots$Area_km2)
    hotspots$Percent_hotspot_sites <- 100 * hotspots$N_sites / sum(hotspots$N_sites)
    hotspots$Percent_hotspot_samples <- 100 * hotspots$N_samples / sum(hotspots$N_samples)
    hotspots$Percent_hotspot_sites[0 == hotspots$Percent_hotspot_sites] <- NA
    hotspots$Percent_hotspot_samples[0 == hotspots$Percent_hotspot_samples] <- NA

    # Log10 scale
    x.range <- range(log10(c(hotspots$Percent_hotspot_area,
                             hotspots$Percent_hotspot_sites,
                             hotspots$Percent_hotspot_samples)), na.rm=TRUE)
    x.offset <- floor(x.range[1])
    xlim <- c(0, (x.range[2] - x.offset) * 1.05)
    col <- rep('grey', nrow(hotspots))
    col[is.na(hotspots$Percent_hotspot_sites)] <- NA
    y <- HorizontalBarplot(log10(hotspots$Percent_hotspot_area) - x.offset,
        xlab=~'% of global hotspot area (bars), of N sites (+) and of N samples ('%*%')',
        show.labels=FALSE, names.arg=NA, xlim=xlim, xaxt='n', col=col, ...)
    x.at <- ceiling(xlim[1]):floor(xlim[2])
    axis(side=1, at=x.at, labels=10^(x.at +x.offset), ...)

    points(log10(hotspots$Percent_hotspot_sites) - x.offset, y, pch=4, ...)    # Crosses
    points(log10(hotspots$Percent_hotspot_samples) - x.offset, y, pch=3, ...)  # Pluses

    # A matrix containing the upper and lower y positions of each predominant realm
    bands <- tapply(y, hotspots[,'Predominant_realm'], range)
    bands <- do.call('rbind', bands)
    bands <- bands[order(bands[,1]),]

    # Plot outside of the normal plotting area
    par(xpd=NA)

    # X axis line across the width of the plot
    segments(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[3], ...)

    # Bands that group predominant realms
    data(realms)
    rect(grconvertX(0, 'ndc', 'user'),
         bands[,1]-diff(y)[1]/2,
         -diff(xlim) * 0.01,
         bands[,2]+diff(y)[1]/2,
         col=as.character(realms$Colour[match(rownames(bands), realms$Realm)]),
         border=NA)

    # Names of predominant realms
    text(grconvertX(0, 'ndc', 'user'),
         apply(bands, 1, function(row) row[1] + diff(row)/2), rownames(bands),
         pos=4, ...)

    # Names of hotspots
    text(0, y, hotspots$Hotspot, pos=2, cex=0.8, ...)
}

PlotHotspotsWildernessAreasScatterCoverage <- function(
  sites, diversity, mar=c(3.5,4.3,0.8,0.5), mgp=c(2.2,0.6,0), las=1, ...) {
    op <- par(mar=mar, mgp=mgp, las=las)
    on.exit(par(op))

    hotspots <- HotspotsTable(sites, diversity)
    hotspots <- hotspots['' != hotspots$Hotspot,]
    hotspots$Name <- hotspots$Hotspot

    # Hotspots colours by predominant realm
    data(realms)
    hotspots$pch <- 19
    hotspots$cex <- 4
    hotspots$Type <- 'H'

    # Wilderness areas
    wilderness.areas <- WildernessAreasTable(sites, diversity)
    wilderness.areas <- wilderness.areas['' != wilderness.areas$Wilderness_area,]
    wilderness.areas$Name <- wilderness.areas$Wilderness_area
    wilderness.areas$pch <- 15
    wilderness.areas$cex <- 4
    wilderness.areas$Type <- 'W'

    # Combine the two tables
    columns <- c('Name', 'Area_km2', 'N_sites', 'Type', 'Predominant_realm',
                 'pch', 'cex')
    to.plot <- rbind(hotspots[,columns], wilderness.areas[,columns])

    to.plot$col <- as.character(
        realms$Alt_colour[match(to.plot$Predominant_realm, realms$Realm)]
    )
    to.plot$col <- paste0(to.plot$col, 'cc')

    unrepresented <- to.plot[0 == to.plot$N_sites,]
    to.plot <- to.plot[!is.na(to.plot$Area_km2) & to.plot$N_sites > 0,]

    # Plot
    to.plot$ID <- 1:nrow(to.plot)
    x <- log10(to.plot$Area_km2)
    y <- log10(to.plot$N_sites)
    plot(x, y, type='n', xlab=~log[10](area)~(km^2), ylab=~log[10](N~Sites))
    points(x, y, col=to.plot$col, pch=to.plot$pch, cex=to.plot$cex)
    text(x, y, to.plot$ID, cex=1)

    # Legend
    # Just those realms that we actually plotted
    realms <- droplevels(realms[realms$Realm %in% to.plot$Predominant_realm,])
    leg <- data.frame(
          text=as.character(realms$Realm),
          col=as.character(realms$Alt_colour),
          pch=rep(19, length(realms$Realm)))
    # Show thick lines in the legend
    legend('bottomright', legend=leg$text, lwd=10, col=as.character(leg$col),
           cex=1.4)

    # Text for caption
    hotspots.caption <- with(to.plot['H' == to.plot$Type,], {
        paste(paste0(ID, ' ', Name), collapse=', ')
    })
    wa.caption <- with(to.plot['W' == to.plot$Type,], {
        paste(paste0(ID, ' ', Name), collapse=', ')
    })
    hotspots.unrepresented <- with(unrepresented['H' == unrepresented$Type,], {
        paste(Name, collapse=', ')
    })
    wa.unrepresented <- with(unrepresented['W' == unrepresented$Type,], {
        paste(Name, collapse=', ')
    })

    invisible(paste0('Hotspots are ', hotspots.caption, ' and ',
                     'HBWAs are ', wa.caption, '. ',
                     'Unrepresented are the hotspots ', hotspots.unrepresented,
                     ' and the HBWA ', wa.unrepresented, '.'))
}
