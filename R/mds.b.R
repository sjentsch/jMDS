mdsClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "mdsClass",
    inherit = mdsBase,
    #### Active bindings ----
    active = list(
        getDta = function() {
            if (is.null(private$.crrDta))
                private$.crrDta <- private$.clnDta()

            return(private$.crrDta)
        },
        getMDS = function() {
            if (is.null(private$.crrMDS))
                private$.crrMDS <- private$.getMDS()

            return(private$.crrMDS)
        }
    ),
    private = list(
        #### Member variables ----
        .crrDta = NULL,
        .crrMDS = NULL,

        #### Init + run functions ----
        .init = function() {
            # initialize table and help / information
            private$.prpCfg()
            private$.genInf()
            private$.addInf()
        },
        .run = function() {
            crrMDS <- self$getMDS
            if (!is.null(crrMDS)) {
                private$.fllCfg()
                for (nmeFig in c("figCfg", "figHst", "figRes", "figShp", "figStr", "figWgh")) self$results[[nmeFig]]$setState(crrMDS)
            }
            private$.genInf()
        },
        # check (and if necessary convert) input data
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
            # the variables to be used in the MDS are converted to numeric
            for (crrVar in mdsVar)
                crrDta[[crrVar]] <- jmvcore::toNumeric(self$data[[crrVar]])
            attr(crrDta, 'row.names') <- seq_along(crrDta[[1]])
            attr(crrDta, 'class') <- 'data.frame'

            return(crrDta)
        },
        # check (and if necessary - i.e., if variables, mode, etc. hav changed - regenerate) MDS estimate
        .getMDS = function() {
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
            } else if (self$options$mdeMDS == "Raw" &&
                       (self$options$dirRaw == "col" && length(self$options$varRaw) >= self$options$dimRaw + 1) ||
                       (self$options$dirRaw == "row" && nrow(crrDta)                >= self$options$dimRaw + 1)) {
                crrMDS <- mdsRaw(crrDta = crrDta, varRaw = self$options$varRaw, nmeRaw = self$options$nmeRaw, xfmRaw = self$options$xfmRaw,
                                 dirRaw = self$options$dirRaw, dimRaw = self$options$dimRaw, lvlRaw = self$options$lvlRaw)
            # [3] input data for determining individual differences (i.e., a series
            #     of symmetric matrices, number of vars * number of individuals)
            } else if (self$options$mdeMDS == "Ind" &&
                       (length(c(self$options$varInd, self$options$nmeInd)) >= self$options$dimInd + 1) &&
                       (nrow(crrDta) %% length(self$options$varInd) == 0) &&
                       (nrow(crrDta) /  length(self$options$varInd) >= 2)) {
                crrMDS <- mdsInd(crrDta = crrDta, varInd = self$options$varInd, nmeInd = self$options$nmeInd,
                                 xfmInd = xfmSnI(self$options), dimInd = self$options$dimInd, lvlInd = self$options$lvlInd)
            } else {
                return(NULL)
            }

            # store estimated solution in tempdir() and return crrMDS
            saveRDS(crrMDS, crrFle)
            
            crrMDS
        },
        # handling general and additional information messages
        .genInf = function() {
            crrMde <- self$options$mdeMDS
            crrInf <- crtInf(self$getMDS, crrMde, self$options[[paste0("lvl", crrMde)]], self$options$dirRaw)
            self$results$genInf$setContent(crrInf)
        },
        .addInf = function() {
             crrInf <- self$results$addInf
             crrMde <- self$options$mdeMDS
             crrInf$setContent(paste0("<p>", paste0(c(addMDS, getVar(paste0("add", crrMde))), collapse = "</p><p>"), "</p>"))
        },
        # prepare and fill configuration table
        .prpCfg = function() {
            crrMde <- self$options$mdeMDS
            crrDim <- self$options[[paste0("dim", crrMde)]]
            # for both symmetric matrices and individual differences, the length of variables plus the name of an eventual name
            # variable determines the number of rows (the name variable is only necessary for sparse matrices, where the number
            # of variables is reduced by 1 - which is compensated for by adding the name variable - length 1)
            # for raw data matrices is becomes slightly more complicated: 
            crrRow <- ifelse(crrMde %in% c("Sym", "Ind"),
                             length(c(self$options[[paste0("var", crrMde)]], self$options[[paste0("nme", crrMde)]])),
                             ifelse(self$options$dirRaw == "col", length(self$options$varRaw), nrow(self$readDataset())))
            if (crrRow < crrDim + 1) return(invisible(NULL))

            crrTbl <- self$results$tblCfg
            nmeClm <- c(sprintf("Dimension %d", seq(crrDim)), rep("SPP", self$options$clmSPP))
            nllRow <- stats::setNames(as.list(rep("", length(nmeClm))), nmeClm)
            for (i in seq_along(nmeClm)) crrTbl$addColumn(name = gsub("Dimension ", "D", nmeClm[i]), title = nmeClm[i], type = 'number')
            for (i in seq(1, crrRow))    crrTbl$addRow(rowKey = i, values = nllRow)
        },
        .fllCfg = function() {
            crrMDS <- self$getMDS
            if (is.null(crrMDS)) return(invisible(NULL))

            crrCfg <- nmeCfg(crrMDS, self$options$dirRaw)
            crrDta <- cbind(data.frame(nmeObj = row.names(crrMDS[[crrCfg]])),
                            as.data.frame(crrMDS[[crrCfg]]),
                            data.frame(SPP = crrMDS[[nmeSPP(crrCfg)]]))
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

            crrCfg <- nmeCfg(crrMDS, self$options$dirRaw)
            crrDta <- cbind(as.data.frame(crrMDS[[crrCfg]]), data.frame(pntSze = bblPnt(crrMDS[[nmeSPP(crrCfg)]], self$options$cfgBbl)))
            crrRng <- c(rndMnM(min(vapply(crrDta[, c("D1", "D2")], min, numeric(1))), 3),
                        rndMnM(max(vapply(crrDta[, c("D1", "D2")], max, numeric(1))), 3))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = D1, y = D2, label = rownames(crrDta))) +
                      ggplot2::geom_point(size = crrDta$pntSze, colour = theme$color[2]) +
                      ggplot2::geom_text(size = 4, vjust = -0.8) +
                      ggplot2::coord_fixed() +
                      ggplot2::labs(x = "Dimension 1", y = "Dimension 2")
            crrFig <- crrFig + if (self$options$cfgInX) ggplot2::scale_x_reverse(limits = rev(crrRng)) else ggplot2::scale_x_continuous(limits = crrRng)
            crrFig <- crrFig + if (self$options$cfgInY) ggplot2::scale_y_reverse(limits = rev(crrRng)) else ggplot2::scale_y_continuous(limits = crrRng)

            print(crrFig + ggtheme)
            TRUE
        },
        # creates a weighted histogram of the dissimilarities, cf. plot(euro.mds, plot.type = "histogram", main = NULL)
        .pltHst = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)

            if (is(crrMDS, "smacofR")) {
                print(pltErr(sprintf("Histogram can not be produced\nfor %s.", crrMDS$model)))
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
                      ggplot2::labs(x = "Dissimiliarity", y = "Frequency") +
                      ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent)
            
            print(crrFig + ggtheme)
            TRUE
        },
        # creates a residual plot, cf. plot(euro.mds, plot.type = "resplot", main = NULL)
        .pltRes = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            
            if (is(crrMDS, "smacofR")) {
                print(pltErr(sprintf("Residual plot can not be produced\nfor %s.", crrMDS$model)))
                return(TRUE)
            }

            if        (is(crrMDS, "smacofB")) {
                crrDta <- data.frame(x = as.numeric(crrMDS$dhat),
                                     y = as.numeric(crrMDS$confdist))
            } else if (is(crrMDS, "smacofID")) {
                crrDta <- data.frame(x = as.numeric(smacof:::sumList(crrMDS$dhat)), y = as.numeric(smacof:::sumList(crrMDS$confdist)))
            }
            crrMax <- rndMnM(max(vapply(crrDta, max, numeric(1))))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x, y = y)) +
                      ggplot2::geom_point(size = 2, color = theme$color[2]) +
                      ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = FALSE, color = theme$color[1]) +
                      ggplot2::coord_fixed() +
                      ggplot2::labs(x = "Normalized Dissimiliarities (d-hats)", y = "Configuration Distances") +
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
                      ggplot2::labs(x = "Observed Dissimilarity", y = "Configuration Distances")

            print(crrFig + ggtheme)
            TRUE
        },
        # create Stress diagram, cf. plot(euro.mds, plot.type="stressplot", main = NULL)
        .pltStr = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)
            
            crrSPP <- crrMDS[[nmeSPP(nmeCfg(crrMDS, self$options$dirRaw))]]
            crrDta <- data.frame(x = seq_along(crrSPP), y = sort(crrSPP, decreasing = TRUE))
            crrFig <- ggplot2::ggplot(crrDta, ggplot2::aes(x = x, y = y, label = rownames(crrDta))) +
                      ggplot2::geom_segment(ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = theme$color[1], linetype = "dotted") +
                      ggplot2::geom_point(size = 1, color = theme$color[2]) +
                      ggplot2::geom_text(size = 4, hjust = 0, vjust = -0.5, angle = 45) +
                      ggplot2::labs(x = "Objects", y = "Stress Proportion (%)") +
                      ggplot2::scale_x_continuous(breaks = seq_along(crrSPP), limits = c(1, max(crrDta$x) + 1), labels = NULL) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(crrSPP) * 1.1))

            print(crrFig + ggtheme)
            TRUE
        },
        .pltWgh = function(image, ggtheme, theme, ...) {
            crrMDS <- image$state
            if (is.null(crrMDS)) return(FALSE)

            if (is(crrMDS, "smacofB") || is(crrMDS, "smacofR")) {
                print(pltErr(sprintf("A weights diagram can not be produced\nfor %s.", crrMDS$model)))
                return(TRUE)
            }
            
            crrDta <- as.data.frame(t(vapply(crrMDS$cweights, diag, numeric(2))), row.names = lstSbj(length(crrMDS$cweights)))
            crrFig <- ggplot2::ggplot(crrDta) +
                      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = D1, y = 0, yend = D2), arrow = ggplot2::arrow()) +
                      ggplot2::geom_text(ggplot2::aes(x = D1, y = D2, label = rownames(crrDta)), hjust = 0, nudge_x = 0.05) +
                      ggplot2::coord_fixed() +
                      ggplot2::labs(x = "Dimension 1", y = "Dimension 2") +
                      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, rndMnM(max(crrDta)))) +
                      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, rndMnM(max(crrDta))))

            print(crrFig + ggtheme)
            TRUE
        }
    )
)
