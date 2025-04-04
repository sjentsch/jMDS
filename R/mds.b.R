#' @importFrom jmvcore .
mdsClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "mdsClass",
    inherit = mdsBase,
    #### Active bindings ----
    active = list(
        getBPD = function() {
            if (is.null(private$.crrBPD))
                private$.crrBPD <- private$.clnBPD()

            return(private$.crrBPD)
        },
        getDta = function() {
            if (is.null(private$.crrDta))
                private$.crrDta <- private$.clnDta()

            return(private$.crrDta)
        },
        getMDS = function() {
            if (is.null(private$.crrMDS))
                private$.crrMDS <- private$.clcMDS()

            return(private$.crrMDS)
        }
    ),
    private = list(
        #### Member variables ----
        .crrBPD = NULL,
        .crrDta = NULL,
        .crrMDS = NULL,

        #### Init + run functions ----
        .init = function() {
            # initialize table and help / information
            private$.prpCfg()
            private$.prpOut()
        },
        .run = function() {
            crrMDS <- self$getMDS
            private$.shwInf()
            if (!is.null(crrMDS)) {
                private$.fllCfg()
                private$.fllOut()
                for (nmeFig in c("figCfg", "figHst", "figRes", "figShp", "figStr", "figWgh")) self$results[[nmeFig]]$setState(crrMDS)
            }
        },
        # check and convert data for biplot
        .clnBPD = function () {
            crrMde <- self$options$mdeMDS
            crrBPD <- self$readDataset()[, self$options[[paste0("bpl", crrMde)]], drop = FALSE]
            if (is.null(self$getMDS) || any(dim(crrBPD) < 1)) return(NULL)
            
            return(as.data.frame(lapply(crrBPD, as.numeric)))
        },
        # check and convert input data
        .clnDta = function() {
            crrMde <- self$options$mdeMDS
            mdsVar <- self$options[[paste0("var", crrMde)]]
            if (length(mdsVar) < self$options[[paste0("dim", crrMde)]] || nrow(self$data) < 1) return(NULL)
            # ensure that the variables have the correct order (as in the dataset; variables might have been
            # added to varSym / varInd in another order and then, the dataset may not be symmetric anymore)
            if (crrMde %in% c("Sym", "Ind")) mdsVar <- intersect(names(self$data), mdsVar)

            crrDta <- list()
            # the name variable is converted to character
            crrVar <- self$options[[paste0("nme", crrMde)]]
            if (length(crrVar) > 0) crrDta[[crrVar]] <- as.character(self$data[[crrVar]])
            # for mode "Ind", there might be an id_Ind variable
            crrVar <- rep(self$options[[paste0("id_", crrMde)]], crrMde == "Ind")
            if (length(crrVar) > 0) crrDta[[crrVar]] <- as.character(self$data[[crrVar]])
            # the variables to be used in the MDS are converted to numeric
            for (crrVar in mdsVar)
                crrDta[[crrVar]] <- jmvcore::toNumeric(self$data[[crrVar]])
            attr(crrDta, 'row.names') <- rownames(self$data)
            attr(crrDta, 'class') <- 'data.frame'

            return(crrDta)
        },
        # check (and if necessary - i.e., if variables, mode, etc. hav changed - regenerate) MDS estimate
        .clcMDS = function() {
            crrDta <- self$getDta
            if (is.null(crrDta)) return(NULL)

            # try to load crrMDS from tempdir() and check it
            crrFle <- file.path(tempdir(), "crrMDS")
            if (file.exists(crrFle)) {
                 crrMDS <- readRDS(crrFle)
                 if (chkMDS(crrMDS, self$options, crrDta)) return(crrMDS)
            }
           
            # if we cannot use the stored / loaded solution (e.g., because settings have changed)
            # we need to calculate it again
            # [1] symmetric input data (e.g., distances, correlations, etc.)
            if        (self$options$mdeMDS == "Sym" &&
                       (length(c(self$options$varSym, self$options$nmeSym)) >= self$options$dimSym + 1) &&
                       (nrow(crrDta) == length(self$options$varSym))) {
                crrMDS <- mdsSym(crrDta = crrDta, varSym = self$options$varSym, nmeSym = self$options$nmeSym,
                                 xfmSym = xfmSnI(self$options), dimSym = self$options$dimSym, lvlSym = self$options$lvlSym)
            # [2] input data as “usual” data matrix (i.e., vars in columns, units / subjects in rows)
            } else if (self$options$mdeMDS == "Rct" &&
                       (self$options$dirRct == "col" && length(self$options$varRct) >= self$options$dimRct + 1) ||
                       (self$options$dirRct == "row" && nrow(crrDta)                >= self$options$dimRct + 1)) {
                crrMDS <- mdsRct(crrDta = crrDta, varRct = self$options$varRct, nmeRct = self$options$nmeRct, xfmRct = self$options$xfmRct,
                                 dirRct = self$options$dirRct, dimRct = self$options$dimRct, lvlRct = self$options$lvlRct)
            # [3] input data for determining individual differences (i.e., a series
            #     of symmetric matrices, number of vars * number of individuals)
            } else if (self$options$mdeMDS == "Ind" &&
                       (length(c(self$options$varInd, self$options$nmeInd)) >= self$options$dimInd + 1) &&
                       (nrow(crrDta) %% length(self$options$varInd) == 0) &&
                       (nrow(crrDta) /  length(self$options$varInd) >= 2)) {
                crrMDS <- mdsInd(crrDta = crrDta, varInd = self$options$varInd, nmeInd = self$options$nmeInd, id_Ind = self$options$id_Ind,
                                 xfmInd = xfmSnI(self$options), dimInd = self$options$dimInd, lvlInd = self$options$lvlInd)
            } else {
                return(NULL)
            }

            # store estimated solution in tempdir() and return crrMDS
            saveRDS(crrMDS, crrFle)
            
            crrMDS
        },
        # create model information and determine which information to show or to hide
        .shwInf = function() {
            crrMde <- self$options$mdeMDS
            crrMDS <- self$getMDS
            if (is.null(crrMDS)) {
                self$results[[paste0("gen", crrMde)]]$setVisible(self$options$shwInf)
                self$results$mdlInf$setVisible(FALSE)
            } else {
                self$results[[paste0("gen", crrMde)]]$setVisible(FALSE)
                crrInf <- c(jmvcore::format(.("Estimated <strong>{m}</strong> (of type \"{t}\") with {o} objects in {i} iterations."),
                                            m = crrMDS$model,
                                            t = getTyp(self$options[[paste0("lvl", crrMde)]]),
                                            o = ifelse(is(crrMDS, "smacofR") && self$options$dirRct == "row", crrMDS$nind, crrMDS$nobj),
                                            i = crrMDS$niter),
                            jmvcore::format(.("Stress value: <strong>{v}</strong>"), v = round(crrMDS$stress, 3)),
                            private$.dcdXfm())
                self$results$mdlInf$setContent(paste(crrInf, collapse = "</p><p>"))
                self$results$mdlInf$setVisible(TRUE)
            }
        },
        # prepare and fill configuration table
        .prpCfg = function() {
            crrMde <- self$options$mdeMDS
            crrDim <- self$options[[paste0("dim", crrMde)]]
            # for both symmetric matrices and individual differences, the length of variables plus the name of an eventual name
            # variable determines the number of rows (the name variable is only necessary for sparse matrices, where the number
            # of variables is reduced by 1 - which is compensated for by adding the name variable - length 1)
            # for data matrices (rectagular) is becomes slightly more complicated: 
            crrRow <- ifelse(crrMde %in% c("Sym", "Ind"),
                             length(c(self$options[[paste0("var", crrMde)]], self$options[[paste0("nme", crrMde)]])),
                             ifelse(self$options$dirRct == "col", length(self$options$varRct), nrow(self$readDataset())))
            if (crrRow < crrDim + 1) return(invisible(NULL))

            crrTbl <- self$results$tblCfg
            nmeClm <- c(jmvcore::format(.("Dimension {d}"), d = seq(crrDim)), rep(.("SPP"), self$options$clmSPP))
            nllRow <- stats::setNames(as.list(rep("", length(nmeClm))), nmeClm)
            for (i in seq_along(nmeClm)) crrTbl$addColumn(name = ifelse(i <= crrDim, sprintf("D%d", i), "SPP"), title = nmeClm[i], type = 'number')
            for (i in seq(1, crrRow))    crrTbl$addRow(rowKey = i, values = nllRow)
        },
        .fllCfg = function() {
            crrMDS <- self$getMDS
            if (is.null(crrMDS)) return(invisible(NULL))

            crrDta <- dtaCfg(crrMDS, self$options$dirRct, "T")
            crrTbl <- self$results$tblCfg
            for (i in seq(nrow(crrDta))) {
                crrTbl$setRow(rowNo = i, as.list(crrDta[i, ]))
            }
        },
        # creates a configuration plot (first two dimensions) or a bubble plot (a configuration plot
        # with the point stress contribution)
        # cf. plot(euro.mds, main = NULL) -> type = "confplot" (default)
        # cf. plot(euro.mds, type = "bubbleplot", main = NULL)
        .pltCfg = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            crrDim <- self$options[[paste0("dim", self$options$mdeMDS)]]
            nteDim <- if (crrDim > 2) jmvcore::format(.("(showing only the first 2 of {d} dimensions)"), d = crrDim) else NULL

            crrDta <- dtaCfg(crrMDS, self$options$dirRct, "P", self$options$cfgBbl, self$options$cfgB4R, theme$color)
            crrRng <- c(rndMnM(min(vapply(crrDta[, c("D1", "D2")], min, numeric(1))), 3),
                        rndMnM(max(vapply(crrDta[, c("D1", "D2")], max, numeric(1))), 3))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = D1, y = D2, label = nmeObj)) +
                      ggplot2::geom_point(size = crrDta$pntSze, color = crrDta$pntClr) +
                      ggplot2::geom_text(size = crrDta$txtSze, color = crrDta$txtClr, vjust = -0.8) +
                      ggplot2::coord_fixed() +
                      ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dotted") +
                      ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dotted") + 
                      ggplot2::labs(x = .("Dimension 1"), y = .("Dimension 2"), caption = nteDim)
            crrFig <- crrFig + if (self$options$cfgInX) ggplot2::scale_x_reverse(limits = rev(crrRng)) else ggplot2::scale_x_continuous(limits = crrRng)
            crrFig <- crrFig + if (self$options$cfgInY) ggplot2::scale_y_reverse(limits = rev(crrRng)) else ggplot2::scale_y_continuous(limits = crrRng)
            if (self$options$cfgBPl && !is.null(self$getBPD)) {
                crrBPl <- smacof::biplotmds(crrMDS, self$getBPD)
                cffBPl <- crrBPl[["coefficients"]] / ceiling(max(crrBPl[["coefficients"]]) / max(abs(crrRng) / 1.2))
                ndgY <- function(v) v + (sign(v) * ifelse(self$options$cfgInY, -1, 1) / 10)
                for (varBPl in colnames(cffBPl)) {
                    crrFig <- crrFig +
                                ggplot2::annotate("segment", x = 0, y = 0, xend = cffBPl[1, varBPl], yend = cffBPl[2, varBPl],
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.20, "cm"), type = "closed"),
                                  size = 0.2, color = theme$color[1])
                    crrFig <- crrFig +
                                ggplot2::annotate("text", x = cffBPl[1, varBPl], y = ndgY(cffBPl[2, varBPl]),
                                  label = sprintf("%s\n(R²: %.3f)", varBPl, crrBPl$R2vec[[varBPl]]),
                                  size = 4, color = theme$color[1])
                }
            }

            print(crrFig + ggtheme)
            TRUE
        },
        # creates a weighted histogram of the dissimilarities, cf. plot(euro.mds, plot.type = "histogram", main = NULL)
        .pltHst = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)

            if (is(crrMDS, "smacofR")) {
                print(pltErr(sprintf(.("Histogram can not be produced\nfor %s."), crrMDS$model)))
                return(TRUE)
            }

            if        (is(crrMDS, "smacofB")) {
                crrDta <- data.frame(x = as.numeric(crrMDS$delta))
            } else if (is(crrMDS, "smacofID")) {
                crrDta <- data.frame(x = sumLst(crrMDS$delta) / length(crrMDS$delta))
            }

            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x)) +
                      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count / sum(count))), bins = 10,
                                              breaks = pretty(crrDta$x, 10), colour = theme$color[1], fill = theme$fill[2]) +
                      ggplot2::labs(x = .("Dissimiliarity"), y = .("Frequency")) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent)
            
            print(crrFig + ggtheme)
            TRUE
        },
        # creates a residual plot, cf. plot(euro.mds, plot.type = "resplot", main = NULL)
        .pltRes = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            
            if (is(crrMDS, "smacofR")) {
                print(pltErr(sprintf(.("Residual plot can not be produced\nfor %s."), crrMDS$model)))
                return(TRUE)
            }

            if        (is(crrMDS, "smacofB")) {
                crrDta <- data.frame(x = as.numeric(crrMDS$dhat),
                                     y = as.numeric(crrMDS$confdist))
            } else if (is(crrMDS, "smacofID")) {
                crrDta <- data.frame(x = sumLst(crrMDS$dhat), y = sumLst(crrMDS$confdist))
            }
            crrMax <- rndMnM(max(vapply(crrDta, max, numeric(1))))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x, y = y)) +
                      ggplot2::geom_point(size = 2, color = theme$color[2]) +
                      ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = FALSE, color = theme$color[1]) +
                      ggplot2::coord_fixed() +
                      ggplot2::labs(x = .("Normalized Dissimiliarities (d-hats)"), y = .("Configuration Distances")) +
                      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, crrMax)) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, crrMax))
            
            print(crrFig + ggtheme)
            TRUE
        },
        # create Shepard diagram, cf. plot(euro.mds, plot.type = "Shepard", main = NULL)
        .pltShp = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            
            if        (is(crrMDS, "smacofB")) {
                crrDta <- data.frame(x = as.numeric(crrMDS$delta),   y = as.numeric(crrMDS$confdist), yf = as.numeric(crrMDS$dhat))
            } else if (is(crrMDS, "smacofR")) {
                crrDta <- data.frame(x = as.numeric(crrMDS$obsdiss), y = as.numeric(crrMDS$confdist), yf = as.numeric(crrMDS$dhat))
            } else if (is(crrMDS, "smacofID")) {
                crrDta <- data.frame(x =     sumLst(crrMDS$delta),   y =     sumLst(crrMDS$confdist), yf =     sumLst(crrMDS$dhat))
            }
            crrDta <- crrDta[order(crrDta$x, crrDta$yf), ]

            crrFit <- sprintf('Non-metric fit, R² = %.3f\nLinear fit, R² = %.3f', 1 - crrMDS$stress ^ 2, cor(crrDta$y, crrDta$yf) ^ 2)
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x, y = y)) +
                      ggplot2::geom_point(size = 1, color = theme$color[2]) +
                      ggplot2::geom_step(ggplot2::aes(x = x, y = yf), color = theme$color[1], direction = "vh") +
                      ggplot2::annotate("text", x = max(crrDta$x) * 0.05, y = max(crrDta$y) * 0.95, label = crrFit, hjust = 0) + 
                      ggplot2::labs(x = .("Observed Dissimilarity"), y = .("Configuration Distances"))

            print(crrFig + ggtheme)
            TRUE
        },
        # create Stress diagram, cf. plot(euro.mds, plot.type="stressplot", main = NULL)
        .pltStr = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            
            crrDta <- data.frame(x = seq_along(crrMDS$spp), y = sort(crrMDS$spp, decreasing = TRUE))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x, y = y, label = rownames(crrDta))) +
                      ggplot2::geom_segment(ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = theme$color[1], linetype = "dotted") +
                      ggplot2::geom_point(size = 1, color = theme$color[2]) +
                      ggplot2::geom_text(size = 4, hjust = 0, vjust = -0.5, angle = 45) +
                      ggplot2::labs(x = .("Objects"), y = .("Stress Proportion (%)")) +
                      ggplot2::scale_x_continuous(breaks = seq_along(crrMDS$spp), limits = c(1, max(crrDta$x) + 1), labels = NULL) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(crrMDS$spp) * 1.2))

            print(crrFig + ggtheme)
            TRUE
        },
        .pltWgh = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)

            if (is(crrMDS, "smacofB") || is(crrMDS, "smacofR")) {
                print(pltErr(sprintf(.("A weights diagram can not be produced\nfor %s."), crrMDS$model)))
                return(TRUE)
            }

            crrDim <- self$options[[paste0("dim", self$options$mdeMDS)]]
            nteDim <- if (crrDim > 2) jmvcore::format(.("(showing only the first 2 of {d} dimensions)"), d = crrDim) else NULL

            crrDta <- as.data.frame(t(vapply(crrMDS$cweights, diag, numeric(crrDim))), row.names = names(crrMDS$delta))
            crrFig <- ggplot2::ggplot(crrDta) +
                      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = D1, y = 0, yend = D2), arrow = ggplot2::arrow()) +
                      ggplot2::geom_text(ggplot2::aes(x = D1, y = D2, label = rownames(crrDta)), hjust = 0, nudge_x = 0.05) +
                      ggplot2::coord_fixed() +
                      ggplot2::labs(x = .("Dimension 1"), y = .("Dimension 2"), caption = nteDim) +
                      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, rndMnM(max(crrDta)))) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, rndMnM(max(crrDta))))

            print(crrFig + ggtheme)
            TRUE
        },
        # writing configuration data to the current data frame
        .prpOut = function() {
            crrMde <- self$options$mdeMDS
            crrDim <- self$options[[paste0("dim", crrMde)]]
            if (self$options$ov_Cfg && length(self$options[[paste0("nme", crrMde)]]) == 0) {
                self$results$ov_Cfg$set(seq(crrDim), sprintf(.("MDS_D%d"), seq(crrDim)),
                                        sprintf(.("MDS Configuration (Dimension %d)"), seq(crrDim)), rep('continuous', crrDim))
            }
        },
        .fllOut = function() {
            crrDta <- self$getDta
            if (is.null(crrDta)) return(invisible(NULL))
            crrMDS <- self$getMDS
            if (is.null(crrMDS)) return(invisible(NULL))

            crrMde <- self$options$mdeMDS
            if (self$options$ov_Cfg && self$results$ov_Cfg$isNotFilled() && length(self$options[[paste0("nme", crrMde)]]) == 0) {
                if        (crrMde == "Sym") {
                    outDta <- as.data.frame(crrMDS$conf)
                } else if (crrMde == "Ind") {
                    outDta <- lst2DF(crrMDS$conf, names(crrMDS$conf[[1]]), mtxTri = FALSE)[, -1]
                } else {
                    return(invisible(NULL))
                }
                if (!is.null(outDta) && nrow(crrDta) == nrow(outDta)) {
                    self$results$ov_Cfg$setRowNums(rownames(crrDta))
                    for (i in seq_along(outDta))
                        self$results$ov_Cfg$setValues(index = i, outDta[, i])
                }
            }
        },
        
        # helper functions ------------------------------------------------------------------------
        # decode transformations - NB: has to be in sync with the transformation operations defined in mds.a.yaml
        .dcdXfm = function() {
            crrMde <- self$options[["mdeMDS"]]
            crrXfm <- self$options[[paste0("xfm", crrMde)]]
            if        (crrMde %in% c("Sym", "Ind")) {
                if        (crrXfm == "none") {
                    .("Matrix contained already distances (i.e., no transformation was applied).")
                } else if (crrXfm == "corr") {
                    .("Before calculating the MDS, the correlations in the data matrix were transformed.")
                } else if (crrXfm == "reverse") {
                    .("Before calculating the MDS, the values in the data matrix were subtracted from the range.")
                } else if (crrXfm == "reciprocal") {
                    .("Before calculating the MDS, the reciprocal of each valiues in the data matrix was calculated.")
                } else if (crrXfm == "ranks") {
                    .("Before calculating the MDS, the values in the data matrix were ranked.")
                } else if (crrXfm == "exp") {
                    .("Before calculating the MDS, the exponential of each valiues in the data matrix was calculated.")
                } else if (crrXfm == "Gaussian") {
                    .("Before calculating the MDS, a Gaussian transformation was applied to the values in the data matrix.")
                } else if (crrXfm == "cooccurrence") {
                    .("Before calculating the MDS, co-occurrences were calculated for the values in the data matrix.")
                } else if (crrXfm == "gravity") {
                    .("Before calculating the MDS, a gravity transformation was applied to the values in the data matrix.")
                } else if (crrXfm == "confusion") {
                    .("Before calculating the MDS, confusion proportions were calculated for the values in the data matrix.")
                } else if (crrXfm == "transition") {
                    .("Before calculating the MDS, transition frequencies were calculated for the values in the data matrix.")
                } else if (crrXfm == "membership") {
                    .("Before calculating the MDS, the membership was calculated for the values in the data matrix.")
                } else if (crrXfm == "probability") {
                    .("Before calculating the MDS, a probability transformation was applied to the values in the data matrix.")
                } else if (is.integer(crrXfm)) {
                        jmvcore::format(.("Before calculating the MDS, the values in the data matrix were subtracted from an integer value ({i})."),
                      i = crrXfm)
                } else {
                    jmvcore::reject(jmvcore::format(.("Invalid transformation {xfm}."), xfm = crrXfm))
                }
            } else if (crrMde %in% c("Rct")) {
                crrDir <- gsub("col", "columns", gsub("row", "rows", self$options[["dirRct"]]))
                dscR2S <- .(" (resulting in a symmetric matrix that afterwards was analyzed using <code>smacofSym</code>)")
                if        (crrXfm == "none") {
                    .("Matrix contained already distances (i.e., no transformation was applied).")
                } else if (crrXfm == "reverse") {
                    .("Before calculating the MDS, the values in the data matrix were subtracted from the range.")
                } else if (crrXfm == "rank") {
                    .("Before calculating the MDS, the values in the data matrix were ranked.")
                } else if (crrXfm %in% c("pearson", "kendall", "spearman")) {
                    fmtStr <- .("Before calculating the MDS, {c}-correlations (over {d}) were calculated and then transformed to distances{i}.")
                    jmvcore::format(fmtStr, c = gsub("\\b([A-Za-z])", "\\U\\1", crrXfm, perl = TRUE), d = crrDir, i = dscR2S)
                } else if (grepl("minkowski_[1-4]", crrXfm)) {
                    fmtStr <- .("Before calculating the MDS, {t} distances ({p}over {d}) were calculated{i}.")
                    jmvcore::format(fmtStr,
                                    t = gsub("minkowski", "Minkowski", gsub("minkowski_2", "Euclidean", gsub("minkowski_1", "Manhattan", crrXfm))),
                                    p = ifelse(grepl("minkowski_[3-4]", crrXfm), paste("power = ", gsub("minkowski_", "", crrXfm), "; "), ""),
                                    d = crrDir, i = dscR2S)
                } else if (grepl("z_minkowski_[1-4]", crrXfm)) {
                    fmtStr <- .("Before calculating the MDS, the data were first z-transformed and then {t} distances ({p}over {d}) were calculated{i}.")
                    jmvcore::format(fmtStr,
                                    t = gsub("z_minkowski", "Minkowski", gsub("z_minkowski_2", "Euclidean", gsub("z_minkowski_1", "Manhattan", crrXfm))),
                                    p = ifelse(grepl("z_minkowski_[3-4]", crrXfm), paste("power = ", gsub("minkowski_", "", crrXfm), "; "), ""),
                                    d = crrDir, i = dscR2S)
                } else if (crrXfm == "binary") {
                    fmtStr <- .("Before calculating the MDS, Jaccard distances (over {d}) were calculated{i}.")
                    jmvcore::format(fmtStr, d = crrDir, i = dscR2S)
                } else if (crrXfm == "z_binary") {
                    fmtStr <- .("Before calculating the MDS, the data were first z-transformed and then Jaccard distances (over {d}) were calculated{i}.")
                    jmvcore::format(fmtStr, d = crrDir, i = dscR2S)
                } else {
                    jmvcore::reject(jmvcore::format(.("Invalid transformation {xfm}."), xfm = crrXfm))
                }
            }
        }

    )
)
