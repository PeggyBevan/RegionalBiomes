PlotErrBar <- function (model, data, responseVar, seMultiplier = 1.96, outDir = NULL, 
                        logLink = "n", catEffects = NULL, contEffects = list(), 
                        contEffectsLabels = NULL, otherCatEffects = list(), otherContEffects = character(0), 
                        forPaper = FALSE, align = FALSE, secdAge = FALSE, xtext.srt = 0, 
                        ylim = NA, yaxp = NULL, order = NULL, rescale = NULL, errbar.cols = NULL,
                        errbar.cap = 0.015, errbar.lwd = 1,
                        pt.pch = NULL, errbar.lty = 1, params = list(), add = FALSE, 
                        offset = 0, plotLabels = TRUE, plotLandUses = TRUE, forcePlotLandUses = FALSE, 
                        plotVertLines = TRUE, compressContEffects = FALSE, cex.txt = NULL, 
                        pt.cex = 1, pt.bg = "white", main = NULL, modelNoInteractions = NULL, 
                        type = "percent", matchLandUses = TRUE, returnColours = FALSE) 
{
  if (!is.null(contEffectsLabels)) {
    if (length(contEffects) != length(contEffectsLabels)) {
      stop("Labels for continuous effects must be of the same number as the effects")
    }
  }
  stopifnot(type %in% c("percent", "response"))
  predefined.factors <- list()
  if (secdAge) {
    predefined.factors$UI <- data.frame(UI = c("Primary Vegetation Minimal use", 
                                               "Primary Vegetation Light use", "Primary Vegetation Intense use", 
                                               "Primary Vegetation Significant use", "Natural Minimal use", 
                                               "Natural Light use", "Natural Intense use", "Natural Significant use", 
                                               "Mature Secondary Vegetation Minimal use", "Mature Secondary Vegetation Light use", 
                                               "Mature Secondary Vegetation Significant use", "Intermediate Secondary Vegetation Minimal use", 
                                               "Intermediate Secondary Vegetation Light use", "Intermediate Secondary Vegetation Significant use", 
                                               "Young Secondary Vegetation Minimal use", "Young Secondary Vegetation Light use", 
                                               "Young Secondary Vegetation Significant use", "Plantation forest Minimal use", 
                                               "Plantation forest Light use", "Plantation forest Intense use", 
                                               "Plantation forest Significant use", "Cropland Minimal use", 
                                               "Cropland Light use", "Cropland Intense use", "Cropland Significant use", 
                                               "Pasture Minimal use", "Pasture Light use", "Pasture Intense use", 
                                               "Pasture Significant use", "Farm Minimal use", "Farm Light use", 
                                               "Farm Intense use", "Farm Significant use", "Agriculture Minimal use", 
                                               "Agriculture Light use", "Agriculture Intense use", 
                                               "Agriculture Significant use", "Urban Minimal use", 
                                               "Urban Light use", "Urban Intense use", "Urban Significant use", 
                                               "Human-dominated Minimal use", "Human-dominated Light use", 
                                               "Human-dominated Intense use", "Human-dominated Significant use"), 
                                        col = c(rep("#66A61E", 4), rep("#66A61E", 4), rep("#147659", 
                                                                                          3), rep("#1B9E77", 3), rep("#8ecfbc", 3), rep("#7570B3", 
                                                                                                                                        4), rep("#E6AB02", 4), rep("#D95F02", 12), rep("#E7298A", 
                                                                                                                                                                                       4), rep("#E7298A", 4)))
    predefined.factors$LandUse <- data.frame(LandUse = c("Primary Vegetation", 
                                                         "Natural", "Mature Secondary Vegetation", "Intermediate Secondary Vegetation", 
                                                         "Young Secondary Vegetation", "Plantation forest", 
                                                         "Cropland", "Pasture", "Farm", "Agriculture", "Urban", 
                                                         "Human-dominated"), col = c("#66A61E", "#66A61E", 
                                                                                     "#147659", "#1B9E77", "#8ecfbc", "#7570B3", "#E6AB02", 
                                                                                     "#D95F02", "#D95F02", "#D95F02", "#E7298A", "#E7298A"))
  } else {
    predefined.factors$UI <- data.frame(UI = c("Primary Vegetation Minimal use", 
                                               "Primary Vegetation Light use", "Primary Vegetation Intense use", 
                                               "Primary Vegetation Significant use", "Natural Minimal use", 
                                               "Natural Light use", "Natural Intense use", "Natural Significant use", 
                                               "Secondary Vegetation Minimal use", "Secondary Vegetation Light use", 
                                               "Secondary Vegetation Intense use", "Secondary Vegetation Significant use", 
                                               "Plantation forest Minimal use", "Plantation forest Light use", 
                                               "Plantation forest Intense use", "Plantation forest Significant use", 
                                               "Cropland Minimal use", "Cropland Light use", "Cropland Intense use", 
                                               "Cropland Significant use", "Pasture Minimal use", 
                                               "Pasture Light use", "Pasture Intense use", "Pasture Significant use", 
                                               "Farm Minimal use", "Farm Light use", "Farm Intense use", 
                                               "Farm Significant use", "Agriculture Minimal use", 
                                               "Agriculture Light use", "Agriculture Intense use", 
                                               "Agriculture Significant use", "Urban Minimal use", 
                                               "Urban Light use", "Urban Intense use", "Urban Significant use", 
                                               "Human-dominated Minimal use", "Human-dominated Light use", 
                                               "Human-dominated Intense use", "Human-dominated Significant use"), 
                                        col = c(rep("#66A61E", 4), rep("#66A61E", 4), rep("#1B9E77", 
                                                                                          4), rep("#7570B3", 4), rep("#E6AB02", 4), rep("#D95F02", 
                                                                                                                                        12), rep("#E7298A", 4), rep("#E7298A", 4)))
    predefined.factors$LandUse <- data.frame(LandUse = c("Primary Vegetation", 
                                                         "Natural", "Null Secondary Vegetation", "Secondary Vegetation", 
                                                         "Null Secondary Vegetation", "Plantation forest", 
                                                         "Cropland", "Pasture", "Farm", "Agriculture", "Urban", 
                                                         "Human-dominated"), col = c("#66A61E", "#66A61E", 
                                                                                     "#1B9E77", "#1B9E77", "#1B9E77", "#7570B3", "#E6AB02", 
                                                                                     "#D95F02", "#D95F02", "#D95F02", "#E7298A", "#E7298A"))
  }
  predefined.factors$UseIntensity <- data.frame(UseIntensity = c("Minimal use", 
                                                                 "Light use", "Intense use"), col = c("#000000", "#000000", 
                                                                                                      "#000000"))
  predefined.factors.2 <- lapply(X = predefined.factors, FUN = function(x) return(x[(!grepl("Natural", 
                                                                                            x[, 1])), ]))
  predefined.factors.2 <- lapply(X = predefined.factors.2, 
                                 FUN = function(x) return(x[(!grepl("Human-dominated", 
                                                                    x[, 1])), ]))
  predefined.factors.2 <- lapply(X = predefined.factors.2, 
                                 FUN = function(x) return(x[(!grepl("Farm", x[, 1])), 
                                                            ]))
  predefined.factors.2 <- lapply(X = predefined.factors.2, 
                                 FUN = function(x) return(x[(!grepl("Agriculture", x[, 
                                                                                     1])), ]))
  predefined.factors.2 <- lapply(X = predefined.factors.2, 
                                 FUN = function(x) return(x[(!grepl("Significant", x[, 
                                                                                     1])), ]))
  labels <- character(0)
  coef.labels <- character(0)
  xVals <- integer(0)
  i <- 1
  for (e in catEffects) {
    if (align & (e %in% names(predefined.factors.2))) {
      eval(substitute(names <- paste(predefined.factors.2$x$x), 
                      list(x = e)))
      labels <- c(labels, names)
      for (n in names) {
        xVals <- c(xVals, i)
        i <- i + 1
      }
      coef.labels <- c(coef.labels, paste(e, names, sep = ""))
    } else if (!matchLandUses) {
      if (class(model)=="lme"){
        eval(substitute(names <- levels(model$data[, 
                                                   e]), list(x = e)))
      } else {
        eval(substitute(names <- levels(model.frame(model)[, 
                                                           e]), list(x = e)))
      }
      
      labels <- c(labels, names)
      for (n in names) {
        xVals <- c(xVals, i)
        i <- i + 1
      }
      coef.labels <- c(coef.labels, paste(e, names, sep = ""))
    } else {
      if (class(model)=="lme"){
        eval(substitute(names <- levels(model$data[, 
                                                   e])[order(match(tolower(levels(model$data[, 
                                                                                             e])), tolower(predefined.factors.2$x$x)))], 
                        list(x = e)))
      } else {
        eval(substitute(names <- levels(model.frame(model)[, 
                                                           e])[order(match(tolower(levels(model.frame(model)[, 
                                                                                                             e])), tolower(predefined.factors.2$x$x)))], 
                        list(x = e)))
      }
      
      labels <- c(labels, names)
      if(!is.null(order)){
        for (n in 1:length(order)) {
          xVals <- c(xVals, i)
          i <- i + 1
        }
      } else {
        for (n in 1:length(names)) {
          xVals <- c(xVals, i)
          i <- i + 1
        }
      }
      coef.labels <- c(coef.labels, paste(e, names, sep = ""))
    }
  }
  o <- match(tolower(coef.labels), tolower(names(fixef(model))))
  y <- fixef(model)[o]
  yplus <- y + se.fixef(model)[o] * seMultiplier
  yminus <- y - se.fixef(model)[o] * seMultiplier
  for (e in catEffects) {
    if (class(model)=="lme"){
      ref.name <- paste(e, levels(model$data[, e])[1], 
                        sep = "")
    } else {
      ref.name <- paste(e, levels(model.frame(model)[, e])[1], 
                        sep = "")
    }
    
    o <- match(tolower(ref.name), tolower(coef.labels))
    y[o] <- 0
    yplus[o] <- 0
    yminus[o] <- 0
  }
  if (!is.null(order)) {
    labels <- labels[order]
    coef.labels <- coef.labels[order]
    y <- y[order]
    yplus <- yplus[order]
    yminus <- yminus[order]
  }
  if (!is.null(rescale)) {
    y <- y + rescale
    yplus <- yplus + rescale
    yminus <- yminus + rescale
  }
  intercept <- fixef(model)["(Intercept)"]
  if (logLink == "e") {
    y <- (exp(y + ifelse(type == "percent", 0, intercept))) - 
      ifelse(type=="percent",0,exp(intercept))
    yplus <- (exp(yplus + ifelse(type == "percent", 0, intercept))) - 
      ifelse(type=="percent",0,exp(intercept))
    yminus <- (exp(yminus + ifelse(type == "percent", 0, 
                                   intercept))) - 
      ifelse(type=="percent",0,exp(intercept))
  } else if (logLink == "10") {
    y <- (10^(y + ifelse(type == "percent", 0, intercept))) - 
      ifelse(type=="percent",0,10^(intercept))
    yplus <- (10^(yplus + ifelse(type == "percent", 0, intercept))) - 
      ifelse(type=="percent",0,10^(intercept))
    yminus <- (10^(yminus + ifelse(type == "percent", 0, 
                                   intercept))) - 
      ifelse(type=="percent",0,10^(intercept))
  } else if (logLink == "inv10") {
    y <- (1/(10^(y)))
    yplus <- (1/(10^(yplus + ifelse(type == "percent", 0, 
                                    intercept))))
    yminus <- (1/(10^(yminus + ifelse(type == "percent", 
                                      0, intercept))))
  } else if (logLink == "b") {
    y <- (1/(1 + exp(-(intercept + y))))
    yplus <- (1/(1 + exp(-(intercept + yplus))))
    yminus <- (1/(1 + exp(-(intercept + yminus))))
  } else if (logLink == "n") {
  } else {
    stop("Error: the specified log link is not supported")
  }
  if ((type == "percent") & (logLink != "n")) {
    if (logLink == "b") {
      y <- ((y/(1/(1 + exp(-(intercept))))) * 100) - 100
      yplus <- ((yplus/(1/(1 + exp(-(intercept))))) * 
                  100) - 100
      yminus <- ((yminus/(1/(1 + exp(-(intercept))))) * 
                   100) - 100
    } else {
      y <- y * 100 - 100
      yplus <- yplus * 100 - 100
      yminus <- yminus * 100 - 100
    }
  } else if ((type == "percent") & (logLink == "n")) {
    y <- ((intercept + y)/intercept) * 100 - 100
    yplus <- ((intercept + yplus)/intercept) * 100 - 100
    yminus <- ((intercept + yminus)/intercept) * 100 - 100
  }
  if (!is.null(modelNoInteractions)) {
    model <- modelNoInteractions
  }
  if (length(contEffects) > 0) {
    if (type != "percent") {
      stop("Currently type must be 'percent' when continuous effects are plotted")
    }
    i2 <- max(xVals + 1)
    for (i in 1:length(contEffects)) {
      if (!is.na(contEffects[i])) {
        if (compressContEffects) {
          xVals <- c(xVals, i2, i2 + (1/3), i2 + (2/3))
          i2 <- i2 + 1
        }
        else {
          xVals <- c(xVals, i2, i2 + 1, i2 + 2)
          i2 <- i2 + 3
        }
        if (is.null(contEffectsLabels)) {
          labels <- c(labels, "", names(contEffects)[i], 
                      "")
        }
        else {
          labels <- c(labels, "", contEffectsLabels[i], 
                      "")
        }
      }
      if (is.na(contEffects[i])) {
      }
      else if (contEffects[i] == 0) {
      }
      else if (contEffects[i] == 1) {
        eval(substitute(newdat <- data.frame(c(min(data$x), 
                                               median(data$x), max(data$x))), list(x = names(contEffects)[i])))
      }
      else if (contEffects[i] == 2) {
        eval(substitute(newdat <- data.frame(c(max(data$x), 
                                               median(data$x), min(data$x))), list(x = names(contEffects)[i])))
      }
      else {
        stop("Continuous effect order must be either NA (don't add to plot), 0 (plot nothing), 1 (ascending) or 2 (descending)")
      }
      if (is.na(contEffects[i])) {
      }
      else if (contEffects[i] == 0) {
        y <- c(y, rep(NA, 3))
        yplus <- c(yplus, rep(NA, 3))
        yminus <- c(yminus, rep(NA, 3))
      }
      else {
        names(newdat) <- names(contEffects)[i]
        for (e in names(contEffects)) {
          if (e != names(contEffects)[i]) {
            newdat[, e] <- median(data[, e])
          }
        }
        for (e in 1:length(otherCatEffects)) {
          newdat[, names(otherCatEffects)[e]] <- factor(otherCatEffects[e], 
                                                        levels = levels(model.frame(model)[, names(otherCatEffects)[e]]))
        }
        if (length(otherContEffects) > 0) {
          for (e in 1:length(otherContEffects)) {
            newdat[, otherContEffects[e]] <- median(data[, 
                                                         otherContEffects[e]])
          }
        }
        for (e in catEffects) {
          newdat[, e] <- factor(levels(model.frame(model)[, 
                                                          e])[1], levels = levels(model.frame(model)[, 
                                                                                                     e]))
        }
        newdat[, names(model.frame(model))[1]] <- 0
        mm <- model.matrix(terms(model), newdat)
        pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)), 
                                              mm))
        newdat$y <- mm %*% fixef(model)
        newdat$yplus <- newdat$y + sqrt(pvar1) * seMultiplier
        newdat$yminus <- newdat$y - sqrt(pvar1) * seMultiplier
        if (logLink == "e") {
          temp.y <- c(0, (((exp(newdat$y[2])/exp(newdat$y[1])) * 
                             100) - 100), (((exp(newdat$y[3])/exp(newdat$y[1])) * 
                                              100) - 100))
          temp.yplus <- c(0, (((exp(newdat$yplus[2])/exp(newdat$y[1])) * 
                                 100) - 100), (((exp(newdat$yplus[3])/exp(newdat$y[1])) * 
                                                  100) - 100))
          temp.yminus <- c(0, (((exp(newdat$yminus[2])/exp(newdat$y[1])) * 
                                  100) - 100), (((exp(newdat$yminus[3])/exp(newdat$y[1])) * 
                                                   100) - 100))
        }
        else if (logLink == "10") {
          temp.y <- c(0, (((10^(newdat$y[2])/10^(newdat$y[1])) * 
                             100) - 100), (((10^(newdat$y[3])/10^(newdat$y[1])) * 
                                              100) - 100))
          temp.yplus <- c(0, (((10^(newdat$yplus[2])/10^(newdat$y[1])) * 
                                 100) - 100), (((10^(newdat$yplus[3])/10^(newdat$y[1])) * 
                                                  100) - 100))
          temp.yminus <- c(0, (((10^(newdat$yminus[2])/10^(newdat$y[1])) * 
                                  100) - 100), (((10^(newdat$yminus[3])/10^(newdat$y[1])) * 
                                                   100) - 100))
        }
        else if (logLink == "inv10") {
          temp.y <- c(0, ((((1/(10^(newdat$y[2])))/(1/(10^(newdat$y[1])))) * 
                             100) - 100), ((((1/(10^(newdat$y[3])))/(1/(10^(newdat$y[1])))) * 
                                              100) - 100))
          temp.yplus <- c(0, ((((1/(10^(newdat$yplus[2])))/(1/(10^(newdat$y[1])))) * 
                                 100) - 100), ((((1/(10^(newdat$yplus[3])))/(1/(10^(newdat$y[1])))) * 
                                                  100) - 100))
          temp.yminus <- c(0, ((((1/(10^(newdat$yminus[2])))/(1/(10^(newdat$y[1])))) * 
                                  100) - 100), ((((1/(10^(newdat$yminus[3])))/(1/(10^(newdat$y[1])))) * 
                                                   100) - 100))
        }
        else if (logLink == "b") {
          intercept <- fixef(model)["(Intercept)"]
          y <- (((1/(1 + exp(-(intercept + y))))/(1/(1 + 
                                                       exp(-(intercept))))) * 100) - 100
          yplus <- (((1/(1 + exp(-(intercept + yplus))))/(1/(1 + 
                                                               exp(-(intercept))))) * 100) - 100
          yminus <- (((1/(1 + exp(-(intercept + yminus))))/(1/(1 + 
                                                                 exp(-(intercept))))) * 100) - 100
        }
        else if (logLink == "n") {
          if (type == "percent") 
            newdat$y <- ((newdat$y + intercept)/intercept) * 
              100 - 100
          temp.y <- c(0, (newdat$y[2] - newdat$y[1]), 
                      (newdat$y[3] - newdat$y[1]))
          temp.yplus <- c(0, (newdat$yplus[2] - newdat$y[1]), 
                          (newdat$yplus[3] - newdat$y[1]))
          temp.yminus <- c(0, (newdat$yminus[2] - newdat$y[1]), 
                           (newdat$yminus[3] - newdat$y[1]))
        }
        else {
          stop("Error: the specified log link is not supported")
        }
        y <- c(y, temp.y)
        yplus <- c(yplus, temp.yplus)
        yminus <- c(yminus, temp.yminus)
      }
    }
  }
  predRange <- max(yplus, na.rm = T) - min(yminus, na.rm = T)
  if (all(is.na(ylim))) {
    plotLims <- c(min(yminus, na.rm = T) - 0.35 * predRange, 
                  max(yplus, na.rm = T) + 0.05 * predRange)
  } else {
    plotLims <- ylim
    predRange <- ylim[2] - ylim[1]
  }
  if (!is.null(outDir)) {
    png(paste(outDir, "/LU effects_", responseVar, ".png", 
              sep = ""), width = 22.86/2.54, height = 12.57/2.54, 
        units = "in", res = 1200)
  }
  if (forPaper) {
    par(mar = c(0.2, 3.5, 0.2, 0.2))
    par(cex.lab = 1)
    par(cex.axis = 0.7)
    cex.pt <- 0.5
    par(mgp = c(2, 1, 0))
    cex.txt <- ifelse(is.null(cex.txt), 0.75, cex.txt)
  } else {
    par(mar = c(1, 5, 1, 1))
    par(cex.lab = 1.8)
    par(cex.axis = 1.5)
    cex.pt <- 1
    par(mgp = c(3.5, 1, 0))
    cex.txt <- ifelse(is.null(cex.txt), 0.8, cex.txt)
  }
  par(las = 1)
  for (p in names(params)) {
    eval(parse(text = gsub("x", p, "par(x=params$x)")))
  }
  if (responseVar != "") {
    if (type != "response") {
      ylabel = paste(responseVar, "difference (%)")
    } else {
      ylabel = paste(responseVar)
    }
  } else {
    ylabel <- ""
  }
  allLevels <- unlist(lapply(predefined.factors, function(x) return(x[, 
                                                                      1])))
  allCols <- unlist(lapply(predefined.factors, function(x) return(x[, 
                                                                    2])))
  plot.cols <- paste(allCols[match(tolower(labels), tolower(allLevels))])
  plot.cols[plot.cols == "NA"] <- "black"
  if (!is.null(errbar.cols)) {
    plot.cols <- errbar.cols
  }
  if (!add) {
    plot(-99999, -99999, xlim = c(0.5, max(xVals) + 0.5), 
         ylim = plotLims, xlab = NA, ylab = ylabel, xaxt = "n", 
         bty = "n", main = main, yaxp = yaxp)
  }
  for (i in 1:length(catEffects)) {
    if ((i != length(catEffects)) | (length(which(!is.na(contEffects))) > 
                                     0)) {
      if (plotVertLines) {
        abline(v = max(grep(catEffects[i], coef.labels)) + 
                 0.5, lwd = 1, col = "dark grey")
      }
    }
  }
  if (any(!is.na(contEffects))) {
    if (any(contEffects > 0)) {
      if (length(which(!is.na(contEffects))) > 0) {
        for (i in 1:length(contEffects)) {
          if (i != length(contEffects)) {
            if (!is.na(contEffects[i])) {
              if (plotVertLines) {
                if (is.null(contEffectsLabels)) {
                  abline(v = which(labels == names(contEffects)[i]) + 
                           1.5, lwd = 1, col = "dark grey")
                }
                else {
                  abline(v = which(labels == contEffectsLabels[i]) + 
                           1.5, lwd = 1, col = "dark grey")
                }
              }
            }
          }
        }
      }
    }
  }
  if ((!add) | forcePlotLandUses) {
    if (TRUE %in% grepl("primary", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("primary", tolower(labels))) - 
               0.5, plotLims[1], max(grep("primary", tolower(labels))) + 
               0.5, plotLims[2], col = "#66A61E33", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("primary", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Primary", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("natural", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("natural", tolower(labels))) - 
               0.5, plotLims[1], max(grep("natural", tolower(labels))) + 
               0.5, plotLims[2], col = "#66A61E33", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("natural", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Natural", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (secdAge) {
      if (TRUE %in% grepl("mature secondary", tolower(labels))) {
        if (!forPaper) {
          rect(min(grep("mature secondary", tolower(labels))) - 
                 0.5, plotLims[1], max(grep("mature secondary", 
                                            tolower(labels))) + 0.5, plotLims[2], col = "#14765933", 
               border = NA)
        }
        if (plotLandUses) {
          text(mean(grep("mature secondary", tolower(labels))) - 
                 0.2, plotLims[1] + 0.05 * predRange, "MSV", 
               col = "black", cex = cex.txt, srt = xtext.srt, 
               xpd = TRUE)
        }
      }
      if (TRUE %in% grepl("intermediate secondary", tolower(labels))) {
        if (!forPaper) {
          rect(min(grep("intermediate secondary", tolower(labels))) - 
                 0.5, plotLims[1], max(grep("intermediate secondary", 
                                            tolower(labels))) + 0.5, plotLims[2], col = "#1B9E7733", 
               border = NA)
        }
        if (plotLandUses) {
          text(mean(grep("intermediate secondary", tolower(labels))) - 
                 0.2, plotLims[1] + 0.05 * predRange, "ISV", 
               col = "black", cex = cex.txt, srt = xtext.srt, 
               xpd = TRUE)
        }
      }
      if (TRUE %in% grepl("young secondary", tolower(labels))) {
        if (!forPaper) {
          rect(min(grep("young secondary", tolower(labels))) - 
                 0.5, plotLims[1], max(grep("young secondary", 
                                            tolower(labels))) + 0.5, plotLims[2], col = "#8ecfbc33", 
               border = NA)
        }
        if (plotLandUses) {
          text(mean(grep("young secondary", tolower(labels))) - 
                 0.2, plotLims[1] + 0.05 * predRange, "YSV", 
               col = "black", cex = cex.txt, srt = xtext.srt, 
               xpd = TRUE)
        }
      }
    }
    else {
      if (TRUE %in% grepl("secondary", tolower(labels))) {
        if (!forPaper) {
          rect(min(grep("secondary", tolower(labels))) - 
                 0.5, plotLims[1], max(grep("secondary", 
                                            tolower(labels))) + 0.5, plotLims[2], col = "#1B9E7733", 
               border = NA)
        }
        if (plotLandUses) {
          text(mean(grep("secondary", tolower(labels))) - 
                 0.2, plotLims[1] + 0.05 * predRange, "Secondary", 
               col = "black", cex = cex.txt, srt = xtext.srt, 
               xpd = TRUE)
        }
      }
    }
    if (TRUE %in% grepl("plantation", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("plantation", tolower(labels))) - 
               0.5, plotLims[1], max(grep("plantation", tolower(labels))) + 
               0.5, plotLims[2], col = "#7570B333", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("plantation", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Plantation", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("cropland", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("cropland", tolower(labels))) - 
               0.5, plotLims[1], max(grep("cropland", tolower(labels))) + 
               0.5, plotLims[2], col = "#E6AB0233", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("cropland", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Cropland", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("pasture", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("pasture", tolower(labels))) - 
               0.5, plotLims[1], max(grep("pasture", tolower(labels))) + 
               0.5, plotLims[2], col = "#D95F0233", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("pasture", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Pasture", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("farm", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("farm", tolower(labels))) - 0.5, 
             plotLims[1], max(grep("farm", tolower(labels))) + 
               0.5, plotLims[2], col = "#D95F0233", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("farm", tolower(labels))) - 0.2, 
             plotLims[1] + 0.05 * predRange, "Farm", col = "black", 
             cex = cex.txt, srt = xtext.srt, xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("agriculture", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("agriculture", tolower(labels))) - 
               0.5, plotLims[1], max(grep("agriculture", 
                                          tolower(labels))) + 0.5, plotLims[2], col = "#D95F0233", 
             border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("agriculture", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Agriculture", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("urban", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("urban", tolower(labels))) - 0.5, 
             plotLims[1], max(grep("urban", tolower(labels))) + 
               0.5, plotLims[2], col = "#E7298A33", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("urban", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Urban", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (TRUE %in% grepl("human", tolower(labels))) {
      if (!forPaper) {
        rect(min(grep("human", tolower(labels))) - 0.5, 
             plotLims[1], max(grep("human", tolower(labels))) + 
               0.5, plotLims[2], col = "#E7298A33", border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("human", tolower(labels))) - 
               0.2, plotLims[1] + 0.05 * predRange, "Human", 
             col = "black", cex = cex.txt, srt = xtext.srt, 
             xpd = TRUE)
      }
    }
    if (any(regexpr("Minimal", labels) == 1)) {
      if (!forPaper) {
        rect(min(which(regexpr("Minimal", labels) == 
                         1)) - 0.5, plotLims[1], max(which(regexpr("Minimal", 
                                                                   labels) == 1)) + 0.5, plotLims[2], col = "#33660033", 
             border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("Minimal", labels)) - 0.2, plotLims[1] + 
               0.05 * predRange, "Minimal", col = "black", 
             cex = cex.txt, srt = xtext.srt, xpd = TRUE)
      }
    }
    if (any(regexpr("Light", labels) == 1)) {
      if (!forPaper) {
        rect(min(which(regexpr("Light", labels) == 1)) - 
               0.5, plotLims[1], max(which(regexpr("Light", 
                                                   labels) == 1)) + 0.5, plotLims[2], col = "#00666633", 
             border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("Light", labels)) - 0.2, plotLims[1] + 
               0.05 * predRange, "Light", col = "black", 
             cex = cex.txt, srt = xtext.srt, xpd = TRUE)
      }
    }
    if (any(regexpr("Intense", labels) == 1)) {
      if (!forPaper) {
        rect(min(which(regexpr("Intense", labels) == 
                         1)) - 0.5, plotLims[1], max(which(regexpr("Intense", 
                                                                   labels) == 1)) + 0.5, plotLims[2], col = "#CC009933", 
             border = NA)
      }
      if (plotLandUses) {
        text(mean(grep("Intense", labels)) - 0.2, plotLims[1] + 
               0.05 * predRange, "Intense", col = "black", 
             cex = cex.txt, srt = xtext.srt, xpd = TRUE)
      }
    }
    if ((plotLabels) & (length(which(plot.cols == "black")) > 
                        0)) {
      text(xVals[which(plot.cols == "black")], plotLims[1] + 
             0.05 * predRange, labels[which(plot.cols == 
                                              "black")], cex = cex.txt, srt = xtext.srt, xpd = TRUE)
    }
    if (logLink == "n") {
      abline(h = 0, col = "dark grey")
    }
    else {
      abline(h = 0, col = "dark grey")
    }
  }
  if (plotLabels) {
    if (length(which(!is.na(contEffects)) > 0)) {
      if ((logLink == "10") | (logLink == "e")) {
        text(xVals[(which(labels == "")[1]:(which(labels == 
                                                    "")[1] + length(which(!is.na(contEffects))) * 
                                              3))] + offset, yminus[which(labels == "")[1]:(which(labels == 
                                                                                                    "")[1] + length(which(!is.na(contEffects))) * 
                                                                                              3)] - (predRange/13), rep(c("L", "M", "H"), 
                                                                                                                        length(which(!is.na(contEffects)))), col = "black", 
             cex = cex.txt, srt = 45)
      } else if (logLink == "inv10") {
        text(xVals[(which(labels == "")[1]:(which(labels == 
                                                    "")[1] + length(which(!is.na(contEffects))) * 
                                              3))] + offset, yplus[which(labels == "")[1]:(which(labels == 
                                                                                                   "")[1] + length(which(!is.na(contEffects))) * 
                                                                                             3)] - (predRange/13), rep(c("L", "M", "H"), 
                                                                                                                       length(which(!is.na(contEffects)))), col = "black", 
             cex = cex.txt, srt = 45)
      } else {
        text(xVals[(which(labels == "")[1]:(which(labels == 
                                                    "")[1] + length(which(!is.na(contEffects))) * 
                                              3))] + offset, yminus[which(labels == "")[1]:(which(labels == 
                                                                                                    "")[1] + length(which(!is.na(contEffects))) * 
                                                                                              3)] - (predRange/13), rep(c("L", "M", "H"), 
                                                                                                                        length(which(!is.na(contEffects)))), col = "black", 
             cex = cex.txt, srt = 45)
      }
    }
  }
  if (is.null(pt.pch)) {
    pt.pch <- 16
  }
  if (is.null(yaxp)) {
    errbar(xVals + offset, y, yplus, yminus, col = plot.cols, 
           errbar.col = plot.cols, add = T, pch = pt.pch, cex = pt.cex, 
           lty = errbar.lty,cap = errbar.cap, lwd = errbar.lwd)
  } else {
    errbar(xVals + offset, y, yplus, yminus, col = plot.cols, 
           errbar.col = plot.cols, add = T, pch = pt.pch, cex = pt.cex, 
           lty = errbar.lty,cap = errbar.cap, lwd = errbar.lwd)
  }
  if ((length(pt.bg) == 1) & (length(labels) > 1)) {
    bg <- rep(pt.bg, length(labels))
    bg[(pt.bg == "match")] <- plot.cols[(pt.bg == "match")]
  } else {
    if (length(pt.bg) != length(labels)) {
      stop("Error point backgrounds must be of same length as number of error bars")
    }
    bg <- pt.bg
    bg[(pt.bg == "match")] <- plot.cols[(pt.bg == "match")]
  }
  points(xVals + offset, y, col = plot.cols, bg = bg, pch = pt.pch, 
         cex = pt.cex)
  if (!is.null(outDir)) {
    dev.off()
  }
  if (returnColours) 
    return(plot.cols)
}


PlotErrBarInter<-function(model,data,responseVar,seMultiplier=1.96,
                          logLink="n",catInteraction=character(0),
                          forPaper=FALSE,secdAge=FALSE,
                          xtext.srt=0,ylim=NA,rescale=NULL,
                          errbar.cols=NULL,
                          errbar.cap = 0.015,errbar.lwd = 1,
                          pt.pch=NULL,
                          pt.bg="white",
                          params=list(),add=FALSE,
                          offset=0,
                          plotLabels=TRUE,
                          plotLandUses=TRUE,
                          forcePlotLandUses=FALSE,
                          cex.txt=NULL,pt.cex=1,
                          main=NULL,type="percent",rescaleRefLevel=TRUE){
  
  if(2!=length(catInteraction)){
    stop("Error: this function only currently supports plotting of one two-way interaction")
  }
  
  if(is.null(pt.pch)){
    pt.pch<-c(16,17,18,15,21,24,23,22,25,4,3,8,9,10,11,12,13,14,7)[
      1:length(levels(model.frame(model)[,catInteraction[2]]))]
  }
  
  nLevels2<-length(levels(data[,catInteraction[2]]))
  
  if (nLevels2 > 5){
    offsets <- seq(from=-0.4,to = 0.4,length.out = nLevels2)+offset
  } else {
    offsets<-seq(from=-(nLevels2-1)/10,to=(nLevels2-1)/10,length.out=nLevels2)+offset
  }
  
  plot.cols<-PlotErrBar(model=model,data=data,responseVar=responseVar,seMultiplier=seMultiplier,
             logLink=logLink,catEffects=catInteraction[1],
             forPaper=forPaper,secdAge=secdAge,xtext.srt=xtext.srt,
             ylim=ylim,rescale=rescale,errbar.cols=errbar.cols,
             errbar.cap = errbar.cap,errbar.lwd = errbar.lwd,
             pt.pch=pt.pch[1],pt.bg = pt.bg[1],
             params=params,add=add,offset=offsets[1],plotLabels=plotLabels,
             plotLandUses=plotLandUses,forcePlotLandUses = forcePlotLandUses,
             cex.txt=cex.txt,pt.cex=pt.cex,main=main,type=type,matchLandUses=FALSE,
             returnColours=TRUE)
  
  fixefs<-fixef(model)
  fixef.names<-names(fixef(model))
  
  var.cov<-vcov(model)
  
  intercept<-fixef(model)['(Intercept)']
  mainEffects1<-c(0,fixefs[setdiff(
    grep(catInteraction[1],fixef.names),grep(":",fixef.names))])
  
  if ((length(pt.bg)==1) & (length(levels(data[,catInteraction[2]]))>1)){
    bg <- rep(pt.bg,length(levels(data[,catInteraction[2]])))
    bg[(pt.bg=="match")] <- plot.cols[(pt.bg=="match")]
  } else {
    if (length(pt.bg) != length(levels(data[,catInteraction[2]]))){
      stop("Error point backgrounds must be of same length as number of error bars")
    }
    bg <- pt.bg
    bg[(pt.bg=="match")] <- plot.cols[(pt.bg=="match")]
  }
  
  invisible(mapply(function(level,offset,pch,pbg){
    
    inter.fe<-c(ifelse(rescaleRefLevel,0,fixefs[paste(catInteraction[2],level,sep="")]),
                mainEffects1[-1]+fixefs[intersect(
                  grep(catInteraction[1],fixef.names),
                  intersect(grep(catInteraction[2],fixef.names),
                  grep(level,fixef.names)))]+
                  ifelse(rescaleRefLevel,0,fixefs[paste(catInteraction[2],level,sep="")]))
    
    covs<-numeric()
    for (l in levels(model.frame(model)[,catInteraction[1]])[-1]){
      covs<-c(covs,var.cov[which(row.names(var.cov)==paste(catInteraction[1],l,sep="")),
                           which(names(var.cov[1,])==paste(catInteraction[1],l,":",
                                                       catInteraction[2],level,sep=""))])
    }
    
    ses1<-se.fixef(model)[setdiff(
      grep(catInteraction[1],fixef.names),grep(":",fixef.names))]
    ses2<-se.fixef(model)[intersect(
      grep(catInteraction[1],fixef.names),
      intersect(grep(catInteraction[2],fixef.names),
                grep(level,fixef.names)))]
    se<-c(0,sqrt(ses1^2+ses2^2+2*covs))
    
    y<-inter.fe
    yplus<-inter.fe+se*seMultiplier
    yminus<-inter.fe-se*seMultiplier

    if (logLink=="e"){
      y<-(exp(y+ifelse(type=="percent",0,intercept)))
      yplus<-(exp(yplus+ifelse(type=="percent",0,intercept)))
      yminus<-(exp(yminus+ifelse(type=="percent",0,intercept)))
    } else if (logLink=="10") {
      y<-(10^(y+ifelse(type=="percent",0,intercept)))
      yplus<-(10^(yplus+ifelse(type=="percent",0,intercept)))
      yminus<-(10^(yminus+ifelse(type=="percent",0,intercept)))
    } else if (logLink=="inv10") {
      y<-(1/(10^(y)))
      yplus<-(1/(10^(yplus+ifelse(type=="percent",0,intercept))))
      yminus<-(1/(10^(yminus+ifelse(type=="percent",0,intercept))))
    } else if (logLink=="b"){
      y<-(1/(1+exp(-(intercept+y))))
      yplus<-(1/(1+exp(-(intercept+yplus))))
      yminus<-(1/(1+exp(-(intercept+yminus))))
    } else if (logLink=="n"){
      
    } else {
      stop("Error: the specified log link is not supported")
    }
    
    if((type=="percent") & (logLink!="n")){
      if(logLink=="b"){
        y<-((y/(1/(1+exp(-(intercept)))))*100)-100
        yplus<-((yplus/(1/(1+exp(-(intercept)))))*100)-100
        yminus<-((yminus/(1/(1+exp(-(intercept)))))*100)-100
      } else {
        y<-y*100-100
        yplus<-yplus*100-100
        yminus<-yminus*100-100
      }
      
      
    }
    
    if ((length(pbg)==1) & (length(inter.fe)>1)){
      bg <- rep(pbg,length(inter.fe))
    } else {
      if (length(pbg) != length(inter.fe)){
        stop("Error point backgrounds must be of same length as number of error bars")
      }
      bg <- pbg
    }
    
    bg[(pbg=="match")] <- plot.cols[(pbg=="match")]
    
    errbar(1:length(inter.fe)+offset,y,yplus,yminus,add=TRUE,col=plot.cols,
           errbar.col=plot.cols,pch=pch,cex=pt.cex,
           cap = errbar.cap,lwd=errbar.lwd)
    
    points(1:length(inter.fe)+offset,y,col=plot.cols,pch=pch,cex=pt.cex,
           bg=bg,cap = errbar.cap,lwd = errbar.lwd)
    
    },as.list(levels(data[,catInteraction[2]])[-1]),
    as.list(offsets[-1]),
    as.list(pt.pch[-1]),
    as.list(bg[-1])))
  
  
}

