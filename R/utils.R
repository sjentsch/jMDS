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

# df2Lst -> see after lst2DF

dst2DF <- function(crrDst = NULL, crrDgn = c()) {
    crrMtx <- as.matrix(crrDst)
    if (!is.null(crrDgn)) diag(crrMtx) <- crrDgn

    as.data.frame(crrMtx)
}

# calculate distance measures for rectangular matrices to be used with smacofSym
dstR4S <- function(crrDta = NULL, varRct = c(), nmeRct = c(), xfmRct = "pearson", dirRct = "col") {
    if        (xfmRct %in% c("pearson", "kendall", "spearman")) {
        smacof::sim2diss(cor(xpsRct(crrDta, varRct, nmeRct, dirRct), method = xfmRct, use = "p"), "corr")
    } else if (xfmRct %in% paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2))) {
        stats::dist(t(xpsRct(crrDta, varRct, nmeRct, dirRct, grepl("^z_", xfmRct))),
                    method = gsub("_[0-9]$", "", gsub("^z_", "", xfmRct)), p = mnkPwr(xfmRct))
    }
}

# calculate distance measures for rectangular matrices to be used with smacofRect
dstR4R <- function(crrDta = NULL, varRct = c(), nmeRct = c(), xfmRct = "none") {
    if        (xfmRct == "none") {
        rowNme(as.matrix(crrDta[, varRct]),                        crrDta[, nmeRct])
    } else if (xfmRct == "rank") {
        rowNme(apply(crrDta[, varRct], 2, rank),                   crrDta[, nmeRct])
    } else if (xfmRct == "reverse") {
        rowNme(apply(crrDta[, varRct], 2, function(x) x - min(x)), crrDta[, nmeRct])
    } else {
        jmvcore::reject("xfmRct must be either \"none\", \"rank\" or \"reverse\".")
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

dtaCfg <- function(crrMDS = NULL, dirRct = "col", outTyp = "T", blnBbl = FALSE, blnBth = TRUE, crrClr = NULL) {
    crrDta <- data.frame()
    # run one loop with dirRct (which is not relevant for symm. and ind. diff.) and
    # run two loops only for rectangular matrices, plots as output, and if show both rows and cols is set
    for (crrDir in c(dirRct, rep(setdiff(c("col", "row"), dirRct), outTyp == "P" && is(crrMDS, "smacofR") && blnBth))) {
        crrCfg <- nmeCfg(crrMDS, crrDir)
        crrDta <- rbind(crrDta, cbind(data.frame(nmeObj = row.names(crrMDS[[crrCfg]])),
                                      as.data.frame(crrMDS[[crrCfg]]),
                                      data.frame(SPP    = crrMDS[[nmeSPP(crrCfg)]],
                                                 pntSze = bblPnt(crrMDS[[nmeSPP(crrCfg)]], blnBbl),
                                                 pntClr = ifelse(is(crrMDS, "smacofR") && blnBth, c(crrDta$pntClr, 0)[1] + 1, 2),
                                                 txtSze = ifelse(is(crrMDS, "smacofR") && blnBth, 3, 4),
                                                 txtClr = ifelse(is(crrMDS, "smacofR") && blnBth, c(crrDta$txtClr, 0)[1] + 1, 1))))
    }
    if (!is.null(crrClr)) crrDta[, c("pntClr", "txtClr")] <- lapply(c("pntClr", "txtClr"), function(C) crrClr[crrDta[, C]])

    crrDta[, -c(rep(ncol(crrDta) - seq(0, 3), outTyp == "T"), rep(ncol(crrDta) - 4, outTyp == "P"))]
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

mdsRct <- function(crrDta = NULL, varRct = c(), nmeRct = c(), xfmRct = "none", dirRct = "col", dimRct = 2, lvlRct = "ordinal") {
    crrArg <- c(list(mdeMDS = "Rct"), as.list(environment()), list(crrHsh = digest::digest(crrDta)))
    # (a) either use the data as they are ("none"), or revert or rank them before using them with smacofRect
    if        (xfmRct %in% c("none", "reverse", "rank")) {
        crrMDS <- smacof::smacofRect(dstR4R(crrDta, varRct, nmeRct, xfmRct),        ndim = dimRct, type = getTyp(lvlRct), ties = getTie(lvlRct))
    # (b) calculate correlations (converted to distances) or distance measures (e.g., euclidean) before using them with smacofSym 
    } else if (xfmRct %in% c("pearson", "kendall", "spearman", paste0(c("", "z_"), rep(c(sprintf("minkowski_%d", seq(4)), "binary"), each = 2)))) {
        crrMDS <- smacof::smacofSym(dstR4S(crrDta, varRct, nmeRct, xfmRct, dirRct), ndim = dimRct, type = getTyp(lvlRct), ties = getTie(lvlRct))
    } else {
        jmvcore::reject(sprintf("xfmRct needs to be one of the following methods: \"%s\"", paste(c("none", "reverse", "rank", "pearson",
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

mnkPwr <- function(xfmMDS = "none") {
    na.omit(c(as.numeric(strsplit(gsub("z_", "", xfmMDS), "_")[[1]][2]), 2))[1]
}

nmeCfg <- function(crrMDS = NULL, dirRct = "col") {
    ifelse(is(crrMDS, "smacofB"), "conf", ifelse(is(crrMDS, "smacofR"), paste0("conf.", dirRct), ifelse(is(crrMDS, "smacofID"), "gspace", "")))
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

rplDsc <- function(crrDsc = NULL, rplDsc = "", rplVar = "", bplVar = NULL) {
    if (nzchar(rplDsc))   crrDsc[["description"]] <- gsub("_RPLDSC_", rplDsc, crrDsc[["description"]])
    if (nzchar(rplVar))   crrDsc[["variables"]]   <- gsub("_RPLVAR_", rplVar, crrDsc[["variables"]])
    if (!is.null(bplVar)) crrDsc[["variables"]]   <-                        c(crrDsc[["variables"]], bplVar)

    crrDsc
}

setDgn <- function(crrDta = NULL, valDgn = 1) {
    crrDta <- as.matrix(crrDta)
    diag(crrDta) <- valDgn
    crrDta
}

sumLst <- function(crrLst = NULL) as.numeric(smacof:::sumList(crrLst))

xpsRct <- function(crrDta = NULL, varRct = c(), nmeRct = c(), dirRct = "col", sclDta = FALSE) {
    if (sclDta) crrDta[, varRct] <- apply(crrDta[, varRct], 2, scale)
    if        (dirRct == "col") {
        rowNme(  as.matrix(crrDta[, varRct]), crrDta[, nmeRct])
    } else if (dirRct == "row") {
        t(rowNme(as.matrix(crrDta[, varRct]), crrDta[, nmeRct]))
    } else {
        jmvcore::reject("dirRct must be either \"col\" or \"row\".")
        invisible(NULL)
    }
}

xfmSnI <- function(crrOpt = NULL) {
    crrMde <- crrOpt$mdeMDS
    ifelse(crrOpt[[paste0("xfm", crrMde)]] == "integer", crrOpt[[paste0("xfi", crrMde)]], crrOpt[[paste0("xfm", crrMde)]])
}
