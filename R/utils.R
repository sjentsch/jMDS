# NB: Functions are in alphabetical order

bblPnt <- function(crrSPP = c(), blnBbl = FALSE) if (blnBbl) crrSPP / length(crrSPP) * 8 else rep(1, length(crrSPP))

chkDgn <- function(crrDta = NULL, xfmSym = "none") {
    if (!is.matrix(crrDta)) crrDta <- as.matrix(crrDta)
    
    (xfmSym %in% c("none") && all(diag(crrDta) %in% c(0, NA))) ||
    (xfmSym %in% c("corr") && all(diag(crrDta) %in% c(1, NA))) ||
    !(xfmSym %in% c("none", "corr"))
}

cnvTrF <- function(crrDta = NULL, xfmSym = "none") {
    if        (all(is.na(crrDta[upper.tri(crrDta)]))) {
        crrDta[upper.tri(crrDta)] <- t(crrDta)[upper.tri(crrDta)]
        if (all(is.na(diag(as.matrix(crrDta))))) crrDta[row(crrDta) == col(crrDta)] <- ifelse(xfmSym %in% c("corr"), 1, 0)
    } else if (all(is.na(crrDta[lower.tri(crrDta)]))) {
        crrDta[lower.tri(crrDta)] <- t(crrDta)[lower.tri(crrDta)]
        if (all(is.na(diag(as.matrix(crrDta))))) crrDta[row(crrDta) == col(crrDta)] <- ifelse(xfmSym %in% c("corr"), 1, 0)
    } else {
        jmvcore::reject("Input data are a triangular matrix, but could not be converted")
        return(invisible(NULL))
    }

    crrDta
}

cnvTrS <- function(crrDta = NULL, varSym = c(), nmeSym = c(), xfmSym = "none") {
    if        (all(is.na(crrDta[, varSym][upper.tri(crrDta[, varSym])])) && all(varSym[-1] == as.character(crrDta[, nmeSym][-nrow(crrDta)]))) {
        crrDta <- setNames(rbind(NA, cbind(crrDta[, varSym], NA)), c(varSym, as.character(crrDta[nrow(crrDta), nmeSym])))
        crrDta[upper.tri(crrDta)] <- t(crrDta)[upper.tri(crrDta)]
        crrDta[row(crrDta) == col(crrDta)] <- ifelse(xfmSym %in% c("corr"), 1, 0)
    } else if (all(is.na(crrDta[, varSym][lower.tri(crrDta[, varSym])])) && all(varSym[-nrow(crrDta)] == as.character(crrDta[, nmeSym][-1]))) {
        crrDta <- setNames(rbind(cbind(NA, crrDta[, varSym]), NA), c(as.character(crrDta[1, nmeSym]), varSym))
        crrDta[lower.tri(crrDta)] <- t(crrDta)[lower.tri(crrDta)]
        crrDta[row(crrDta) == col(crrDta)] <- ifelse(xfmSym %in% c("corr"), 1, 0)
    } else {
        jmvcore::reject("Input data are a sparse triangular matrix, but it could not be converted")
        return(invisible(NULL))
    }
    row.names(crrDta) <- names(crrDta)
    
    crrDta
}

# calculate distance measures for raw data to be used with smacofSym
dstR4S <- function(crrDta = NULL, varRaw = c(), nmeRaw = c(), xfmRaw = "pearson", dirRaw = "col") {
    if        (xfmRaw %in% c("pearson", "kendall", "spearman")) {
        smacof::sim2diss(cor(xpsRaw(crrDta, varRaw, nmeRaw, dirRaw), method = xfmRaw, use = "p"), "corr")
    } else if (xfmRaw %in% paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2))) {
        stats::dist(t(xpsRaw(crrDta, varRaw, nmeRaw, dirRaw, grepl("^z_", xfmRaw))),
                    method = gsub("_[0-9]$", "", gsub("^z_", "", xfmRaw)), p = mnkPwr(xfmRaw))
    }
}

# calculate distance measures for raw data to be used with smacofRect
dstR4R <- function(crrDta = NULL, varRaw = c(), nmeRaw = c(), xfmRaw = "none") {
    if        (xfmRaw == "none") {
        rowNme(as.matrix(crrDta[, varRaw]),                        crrDta[, nmeRaw])
    } else if (xfmRaw == "rank") {
        rowNme(apply(crrDta[, varRaw], 2, rank),                   crrDta[, nmeRaw])
    } else if (xfmRaw == "reverse") {
        rowNme(apply(crrDta[, varRaw], 2, function(x) x - min(x)), crrDta[, nmeRaw])
    } else {
        jmvcore::reject("xfmRaw must be either \"none\", \"rank\" or \"reverse\".")
        invisible(NULL)
    }
}

dstSym <- function(crrDta = NULL, varSym = c(), nmeSym = c(), xfmSym = "none") {
    # checks whether the input data frame is symmetric, transform if not
    if (isSymmetric(unname(as.matrix(crrDta)))) {
        crrDta <- as.matrix(crrDta)
        rownames(crrDta) <- varSym
    } else {
       if (diff(dim(crrDta)) == 0 && !is.null(varSym) && length(varSym) == dim(crrDta)[1]) {
           crrDta <- cnvTrF(crrDta, xfmSym)
           rownames(crrDta) <- varSym
       } else if (diff(dim(crrDta)) == 1 && !is.null(nmeSym)) {
           crrDta <- cnvTrS(crrDta, varSym, nmeSym, xfmSym)
       } else {
           jmvcore::reject("Input data are neither a symmetric nor a triangular matrix")
           return(invisible(NULL))
       }
    }

    # convert the data (first to a matrix and then to distances) and return them
    if (!chkDgn(crrDta, xfmSym)) {
        jmvcore::reject(sprintf(paste("Input data are a symmetric matrix (or were converted into one), but values in the",
                                      "main diagonal have to be either %d or NA (indicating that it doesn't contain %s)."),
                                      ifelse(xfmSym %in% c("corr"), 1, 0),
                                      ifelse(xfmSym %in% c("corr"), "correlations", "distances")))
        return(invisible(NULL))
    }

    # convert data.frame as distances (if they are already distances),
    # or calculate distances if the data are similarities
    # afterwards calculate smacofSym
    if (xfmSym == "none") {
        as.dist(crrDta) 
    } else {
        smacof::sim2diss(crrDta, method = xfmSym, to.dist = TRUE)
    }
}

mnkPwr <- function(xfmRaw = "none") {
    na.omit(c(as.numeric(strsplit(gsub("z_", "", xfmRaw), "_")[[1]][2]), 2))[1]
}

getVar <- function(inpStr = "") {
    eval(parse(text = inpStr))
}

lstSbj <- function(numSbj = 1) {
    sprintf(paste0("S_%0", ceiling(log10(numSbj + 1)), "d"), seq(numSbj))
}

mdsInd <- function(crrDta = NULL, varInd = c(), nmeInd = c(), xfmInd = "none", dimInd = 2, lvlInd = "ordinal") {
    lngInd <- length(varInd)
    # check whether the length of the data frame is a multiple of varInd
    if (nrow(crrDta) %% lngInd != 0) {
        jmvcore::reject(sprintf(paste("The number of rows in the input data (%d) needs to be a multiple of the number of",
                                      "variables to be included in the MDS. Please check whether there are additional rows",
                                      "in the dataset (e.g., empty rows at the end)."), nrow(crrDta), lngInd))
        return(invisible(NULL))
    }

    numSbj <- nrow(crrDta) / lngInd
    crrLst <- vector(mode = "list", length = numSbj)
    for (crrSbj in seq(numSbj)) {
        crrLst[[crrSbj]] <- dstSym(crrDta[crrSbj * lngInd - seq(lngInd - 1, 0), ], varInd, nmeInd, xfmInd)
    }

    smacof::smacofIndDiff(crrLst, ndim = dimInd, type = lvlInd)
}

mdsRaw <- function(crrDta = NULL, varRaw = c(), nmeRaw = c(), xfmRaw = "none", dirRaw = "col", dimRaw = 2, lvlRaw = "ordinal") {
    # (a) either use the data as they are ("none"), or revert or rank them before using them with smacofRect
    if        (xfmRaw %in% c("none", "reverse", "rank")) {
        smacof::smacofRect(dstR4R(crrDta, varRaw, nmeRaw, xfmRaw),        ndim = dimRaw, type = lvlRaw)
    # (b) calculate correlations (converted to distances) or distance measures (e.g., euclidean) before using them with smacofSym 
    } else if (xfmRaw %in% c("pearson", "kendall", "spearman", paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2)))) {
        smacof::smacofSym(dstR4S(crrDta, varRaw, nmeRaw, xfmRaw, dirRaw), ndim = dimRaw, type = lvlRaw)
    } else {
        jmvcore::reject(sprintf("xfmRaw needs to be one of the following methods: \"%s\"", paste(c("none", "reverse", "rank", "pearson",
          "kendall", "spearman", paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2))), collapse = "\", \"")))
    }
}

mdsSym <- function(crrDta = NULL, varSym = c(), nmeSym = c(), xfmSym = "none", dimSym = 2, lvlSym = "ordinal") {
    smacof::smacofSym(dstSym(crrDta, varSym, nmeSym, xfmSym), ndim = dimSym, type = lvlSym)
}

nmeCfg <- function(crrMDS = NULL, dirRaw = "col") {
    ifelse(is(crrMDS, "smacofB"), "conf", ifelse(is(crrMDS, "smacofR"), paste0("conf.", dirRaw), ifelse(is(crrMDS, "smacofID"), "gspace", "")))
}

nmeSPP <- function(crrCfg = c()) gsub("conf|gspace", "spp", crrCfg)

pltErr <- function(errMsg = "") {
    ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1, size = 8, label = errMsg) +
    ggplot2::theme_void()
}

rndMnM <- function(crrExt = 0, addMlt = 1) {
    crrMlt <- 10 ^ (floor(log10(abs(crrExt))) - 1)
    if        (crrExt > 0) {
        ceiling(crrExt / crrMlt + addMlt) * crrMlt
    } else if (crrExt < 0) {
        floor(crrExt   / crrMlt - addMlt) * crrMlt
    } else {
        crrExt
    }
}

rowNme <- function(crrDta = NULL, valNme = c()) {
    row.names(crrDta) <- if (length(valNme) > 0) valNme else lstSbj(nrow(crrDta))
    crrDta
}

sumLst <- function(crrLst = NULL) as.numeric(smacof:::sumList(crrLst))

xpsRaw <- function(crrDta = NULL, varRaw = c(), nmeRaw = c(), dirRaw = "col", sclDta = FALSE) {
    if (sclDta) crrDta[, varRaw] <- apply(crrDta[, varRaw], 2, scale)
    if        (dirRaw == "col") {
        rowNme(  as.matrix(crrDta[, varRaw]), crrDta[, nmeRaw])
    } else if (dirRaw == "row") {
        t(rowNme(as.matrix(crrDta[, varRaw]), crrDta[, nmeRaw]))
    } else {
        jmvcore::reject("dirRaw must be either \"col\" or \"row\".")
        invisible(NULL)
    }
}
