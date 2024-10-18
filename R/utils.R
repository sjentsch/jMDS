# NB: Functions are in alphabetical order

attLbl <- function(crrDta = NULL, crrVar = c()) {
    for (C in names(crrDta)) attr(crrDta[[C]], "jmv-desc") <- crrVar[[C]]
    
    crrDta
}

bblPnt <- function(crrSPP = c(), blnBbl = FALSE) if (blnBbl) crrSPP / length(crrSPP) * 8 else rep(1, length(crrSPP))

chkDgn <- function(crrDta = NULL, xfmSym = "none") {
    if (!is.matrix(crrDta)) crrDta <- as.matrix(crrDta)
    
    (xfmSym %in% c("none") && all(diag(crrDta) %in% c(0, NA))) ||
    (xfmSym %in% c("corr") && all(diag(crrDta) %in% c(1, NA))) ||
    !(xfmSym %in% c("none", "corr"))
}

chkMDS <- function(crrMDS = NULL, crrOpt = NULL, crrDta = NULL) {
    crrArg <- attr(crrMDS, "crrArg")
    eqlMDS <- crrArg[["crrHsh"]] == digest::digest(crrDta)
    for (N in setdiff(names(crrArg), "crrHsh")) eqlMDS <- eqlMDS && identical(crrArg[[N]], crrOpt[[N]])
    
    eqlMDS
}

cnvSps <- function(crrDta = NULL, varSym = c(), nmeSym = c(), xfmSym = "none") {
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

cnvTrn <- function(crrDta = NULL, xfmSym = "none") {
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

crtInf <- function(crrMDS = NULL, crrMde = "Sym", crrLvl = "", crrDrR = "col") {
    if (!is.null(crrMDS)) {
        sprintf(paste("<p>Estimated <strong>%s</strong> (of type \"%s\") with %d objects in %d iterations.</p>",
                      "<p>Stress-1 value: <strong>%.4f</strong></p><p>%s</p>"),
                crrMDS$model, getTyp(crrLvl), ifelse(is(crrMDS, "smacofR") && crrDrR == "row", crrMDS$nind, crrMDS$nobj),
                crrMDS$niter, crrMDS$stress, dcdXfm(attr(crrMDS, "crrArg")))
    } else {
        paste0("<p>", paste0(c(genMDS, getVar(paste0("gen", crrMde))), collapse = "</p><p>"), "</p>")
    }
}

# decode transformations - has to be in sync with the
# transformation operations defined in mds.a.yaml
dcdXfm <- function(crrArg = NULL) {
    crrMde <- crrArg[["mdeMDS"]]
    crrXfm <- crrArg[[paste0("xfm", crrMde)]]
    if        (crrMde %in% c("Sym", "Ind")) {
        if        (crrXfm == "none") {
            "Matrix contained already distances (i.e., no transformation was applied)."
        } else if (crrXfm == "corr") {
            "Before calculating the MDS, the correlations in the data matrix were transformed."
        } else if (crrXfm == "reverse") {
            "Before calculating the MDS, the values in the data matrix were subtracted from the range."
        } else if (crrXfm == "reciprocal") {
            "Before calculating the MDS, the reciprocal of each valiues in the data matrix was calculated."
        } else if (crrXfm == "ranks") {
            "Before calculating the MDS, the values in the data matrix were ranked."
        } else if (crrXfm == "exp") {
            "Before calculating the MDS, the exponential of each valiues in the data matrix was calculated."
        } else if (crrXfm == "Gaussian") {
            "Before calculating the MDS, a Gaussian transformation was applied to the values in the data matrix."
        } else if (crrXfm == "cooccurrence") {
            "Before calculating the MDS, co-occurrences were calculated for the values in the data matrix."
        } else if (crrXfm == "gravity") {
            "Before calculating the MDS, a gravity transformation was applied to the values in the data matrix."
        } else if (crrXfm == "confusion") {
            "Before calculating the MDS, confusion proportions were calculated for the values in the data matrix."
        } else if (crrXfm == "transition") {
            "Before calculating the MDS, transition frequencies were calculated for the values in the data matrix."
        } else if (crrXfm == "membership") {
            "Before calculating the MDS, the membership was calculated for the values in the data matrix."
        } else if (crrXfm == "probability") {
            "Before calculating the MDS, a probability transformation was applied to the values in the data matrix."
        } else if (is.integer(crrXfm)) {
            sprintf("Before calculating the MDS, the values in the data matrix were subtracted from an integer value (%d).", crrXfm)
        } else {
            jmvcore::reject(sprintf("Invalid transformation %s.", crrXfm))
            return(c())
        }
    } else if (crrMde %in% c("Raw")) {
        crrDir <- gsub("col", "columns", gsub("row", "rows", crrArg[["dirRaw"]]))
        dscR2S <- "(resulting in a symmetric matrix that afterwards was analyzed using <code>smacofSym</code>)"
        if        (crrXfm == "none") {
            "Matrix contained already distances (i.e., no transformation was applied)."
        } else if (crrXfm == "reverse") {
            "Before calculating the MDS, the values in the data matrix were subtracted from the range."
        } else if (crrXfm == "rank") {
            "Before calculating the MDS, the values in the data matrix were ranked."
        } else if (crrXfm %in% c("pearson", "kendall", "spearman")) {
            sprintf("Before calculating the MDS, %s-correlations (over %s) were calculated and then transformed to distances %s.",
                    gsub("\\b([A-Za-z])", "\\U\\1", crrXfm, perl = TRUE), crrDir, dscR2S)
        } else if (grepl("minkowski_[1-4]", crrXfm)) {
            sprintf("Before calculating the MDS, %s distances (%sover %s) were calculated%s.",
                    gsub("minkowski", "Minkowski", gsub("minkowski_2", "Euclidean", gsub("minkowski_1", "Manhattan", crrXfm))),
                    ifelse(grepl("minkowski_[3-4]", crrXfm), paste("power = ", gsub("minkowski_", "", crrXfm)), ""),
                    crrDir, dscR2S)
        } else if (grepl("z_minkowski_[1-4]", crrXfm)) {
            sprintf("Before calculating the MDS, the data were first z-transformed and then %s distances (%sover %s) were calculated %s.",
                    gsub("z_minkowski", "Minkowski", gsub("z_minkowski_2", "Euclidean", gsub("z_minkowski_1", "Manhattan", crrXfm))),
                    ifelse(grepl("z_minkowski_[3-4]", crrXfm), paste("power = ", gsub("minkowski_", "", crrXfm)), ""),
                    crrDir, dscR2S)
        } else if (crrXfm == "binary") {
            sprintf("Before calculating the MDS, Jaccard distances (over %s) were calculated%s.", crrDir, dscR2S)
        } else if (crrXfm == "z_binary") {
            sprintf("Before calculating the MDS, the data were first z-transformed and then Jaccard distances (over %s) were calculated %s.",
                    crrDir, dscR2S)
        } else {
            jmvcore::reject(sprintf("Invalid transformation %s.", crrXfm))
            return(c())
        }
    }
}

# df2Lst -> see after lst2DF

dst2DF <- function(crrDst = NULL, crrDgn = c()) {
    crrMtx <- as.matrix(crrDst)
    if (!is.null(crrDgn)) diag(crrMtx) <- crrDgn

    as.data.frame(crrMtx)
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
    # remove other columns than varSym and nmeSym
    crrDta <- crrDta[, c(nmeSym, varSym)]
    # checks whether the input data frame is symmetric, transform if not
    if (isSymmetric(unname(as.matrix(crrDta)))) {
        crrDta <- as.matrix(crrDta)
        rownames(crrDta) <- varSym
    } else {
       if (diff(dim(crrDta)) == 0 && !is.null(varSym) && length(varSym) == dim(crrDta)[1]) {
           crrDta <- cnvTrn(crrDta, xfmSym)
           rownames(crrDta) <- varSym
       } else if (diff(dim(crrDta)) == 1 && !is.null(nmeSym)) {
           crrDta <- cnvSps(crrDta, varSym, nmeSym, xfmSym)
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

getID  <- function(crrID = NULL) {
    if (length(crrID) == 0) return(c())
    
    if (!all(diff(table(crrID)) == 0)) {
        jmvcore::reject(sprintf(paste("The column %s (chosen as ID variable) has not the same number of entries for all individuals.",
                                      "This might be due to additional empty lines in the data file."), attr(crrID, "name")))
        return(c())
    }
    if (sum(c(crrID[-1], crrID[1]) != crrID) != length(unique(crrID))) {
        jmvcore::reject(sprintf(paste("The values in the column %s (chosen as ID variable) need to be consecutive for each individual.",
                                      "That is, all values for the first individual need to come first, then all values for the second",
                                      "individual, and so on."), attr(crrID, "name")))
        return(c())
    }

    as.character(unique(crrID))
}

getTie <- function(lvlStr = "") {
    c(na.omit(strsplit(lvlStr, "_")[[1]][2]), "primary")[1]
}

getTyp <- function(lvlStr = "") {
    strsplit(lvlStr, "_")[[1]][1]
}

getVar <- function(inpStr = "") {
    eval(parse(text = inpStr))
}

lstSbj <- function(numSbj = 1) {
    sprintf(paste0("S_%0", ceiling(log10(numSbj + 1)), "d"), seq(numSbj))
}

lst2DF <- function(crrLst = NULL, nmeVar = c(), mtxTri = TRUE, mtxSps = FALSE, valDgn = NA) {
    crrLst <- lapply(crrLst, function(M) { M <- as.matrix(M); if (mtxTri || mtxSps) M[upper.tri(M)] <- NA; as.data.frame(M) })
    if (is.null(names(crrLst)))
        names(crrLst) <- lstSbj(length(crrLst))
    if (!is.null(nmeVar))
        crrLst <- lapply(crrLst, function(M) { names(M) <- nmeVar; M })
    if (!is.na(valDgn))
        crrLst <- lapply(crrLst, function(M) { as.data.frame(setDgn(M, valDgn)) })
    if (mtxSps)
        crrLst <- lapply(crrLst, function(M) { M <- as.matrix(M); diag(M) <- NA; cbind(data.frame(Name = colnames(M)), as.data.frame(M))[-1, seq(nrow(M))] })
    crrLsI <- nrow(crrLst[[1]])
    crrDta <- setNames(as.data.frame(matrix(NA, nrow = length(crrLst) * crrLsI, ncol = ncol(crrLst[[1]]) + 1)), c("ID", names(crrLst[[1]])))
    for (C in seq_along(crrLst))
        crrDta[C * crrLsI - seq(crrLsI - 1, 0), ] <- cbind(data.frame(ID = names(crrLst)[C]), crrLst[[C]])
    
    crrDta
}

df2Lst <- function(crrDta = NULL, crrVar = c(), crrNme = c(), crrID = c(), crrXfm = "none") {
    crrNmV <- length(crrVar)
    crrNmR <- nrow(crrDta) 
    # check whether the length of the data frame is a multiple of crrVar
    if (crrNmR %% crrNmV != 0 && crrNmR / crrNmV >= 2) {
        jmvcore::reject(sprintf(paste("The number of rows in the input data (%d) needs to be a multiple of the number of",
                                      "variables to be included in the MDS. Please check whether there are additional rows",
                                      "in the dataset (e.g., empty rows at the end)."), crrNmR, crrNmV))
        return(invisible(NULL))
    }

    crrNmS <- crrNmR / crrNmV
    # if crrID is empty, assign create a vector with IDs (S_...)
    crrLst <- setNames(vector(mode = "list", length = crrNmS), if (is.null(crrID)) lstSbj(crrNmS) else getID(crrDta[, crrID]))
    for (crrSbj in seq(crrNmS)) {
        crrLst[[crrSbj]] <- dstSym(crrDta[crrSbj * crrNmV - seq(crrNmV - 1, 0), ], crrVar, crrNme, crrXfm)
    }

    crrLst
}

mdsInd <- function(crrDta = NULL, varInd = c(), nmeInd = c(), id_Ind = c(), xfmInd = "none", dimInd = 2, lvlInd = "ordinal") {
    crrArg <- c(list(mdeMDS = "Ind"), as.list(environment()), list(crrHsh = digest::digest(crrDta)))

    # ensure that the variables in varInd appear in the same order as in the data set
    varInd <- intersect(names(crrDta), varInd)

    crrMDS <- smacof::smacofIndDiff(df2Lst(crrDta, crrVar = varInd, crrID = id_Ind, crrXfm = xfmInd),
                                    ndim = dimInd, type = getTyp(lvlInd), ties = getTie(lvlInd))
    attr(crrMDS, "crrArg") <- crrArg[setdiff(names(crrArg), "crrDta")]

    crrMDS
}

mdsRaw <- function(crrDta = NULL, varRaw = c(), nmeRaw = c(), xfmRaw = "none", dirRaw = "col", dimRaw = 2, lvlRaw = "ordinal") {
    crrArg <- c(list(mdeMDS = "Raw"), as.list(environment()), list(crrHsh = digest::digest(crrDta)))
    # (a) either use the data as they are ("none"), or revert or rank them before using them with smacofRect
    if        (xfmRaw %in% c("none", "reverse", "rank")) {
        crrMDS <- smacof::smacofRect(dstR4R(crrDta, varRaw, nmeRaw, xfmRaw),        ndim = dimRaw, type = getTyp(lvlRaw), ties = getTie(lvlRaw))
    # (b) calculate correlations (converted to distances) or distance measures (e.g., euclidean) before using them with smacofSym 
    } else if (xfmRaw %in% c("pearson", "kendall", "spearman", paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2)))) {
        crrMDS <- smacof::smacofSym(dstR4S(crrDta, varRaw, nmeRaw, xfmRaw, dirRaw), ndim = dimRaw, type = getTyp(lvlRaw), ties = getTie(lvlRaw))
    } else {
        jmvcore::reject(sprintf("xfmRaw needs to be one of the following methods: \"%s\"", paste(c("none", "reverse", "rank", "pearson",
          "kendall", "spearman", paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2))), collapse = "\", \"")))
    }
    attr(crrMDS, "crrArg") <- crrArg[setdiff(names(crrArg), "crrDta")]
    
    crrMDS
}

mdsSym <- function(crrDta = NULL, varSym = c(), nmeSym = c(), xfmSym = "none", dimSym = 2, lvlSym = "ordinal") {
    crrArg <- c(list(mdeMDS = "Sym"), as.list(environment()), list(crrHsh = digest::digest(crrDta)))
    # ensure that the variables in varSym appear in the same order as in the data set
    varSym <- intersect(names(crrDta), varSym)
    crrMDS <- smacof::smacofSym(dstSym(crrDta, varSym, nmeSym, xfmSym), ndim = dimSym, type = getTyp(lvlSym), ties = getTie(lvlSym))
    attr(crrMDS, "crrArg") <- crrArg[setdiff(names(crrArg), "crrDta")]
    
    crrMDS   
}

mnkPwr <- function(xfmRaw = "none") {
    na.omit(c(as.numeric(strsplit(gsub("z_", "", xfmRaw), "_")[[1]][2]), 2))[1]
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

rplDsc <- function(crrDsc = NULL, rplDsc = "", rplVar = "") {
    if (nzchar(rplDsc)) crrDsc[["description"]] <- gsub("_RPLDSC_", rplDsc, crrDsc[["description"]])
    if (nzchar(rplVar)) crrDsc[["variables"]]   <- gsub("_RPLVAR_", rplVar, crrDsc[["variables"]])

    crrDsc
}

setDgn <- function(crrDta = NULL, valDgn = 1) {
    crrDta <- as.matrix(crrDta)
    diag(crrDta) <- valDgn
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

xfmSnI <- function(crrOpt = NULL) {
    crrMde <- crrOpt$mdeMDS
    ifelse(crrOpt[[paste0("xfm", crrMde)]] == "integer", crrOpt[[paste0("xfi", crrMde)]], crrOpt[[paste0("xfm", crrMde)]])
}
