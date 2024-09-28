mdsClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "mdsClass",
    inherit = mdsBase,
    private = list(
        .crrDta = NULL,
        .crrMDS = NULL,

        .init = function() {
            private$.crrDta <- NULL
            private$.crrMDS <- NULL

            if (is.null(self$options[[paste0("var", self$options$mdeMDS)]])) return(invisible(NULL))

            # initialize table and help / information
            # [1] symmetric input data (e.g., distances, correlations, etc.)
            if        (self$options$mdeMDS == "Sym") {
                private$.prpCfg(crrDim = self$options$dimSym, numRow = length(c(self$options$varSym, self$options$nmeSym)))
            # [2] input data as “usual” data matrix (i.e., vars in columns, units / subjects in rows)
            } else if (self$options$mdeMDS == "Raw") {
                crrRow <- ifelse(self$options$dirRaw == "col", length(self$options$varRaw), nrow(self$readDataset()))
                private$.prpCfg(crrDim = self$options$dimRaw, numRow = crrRow)
            # [3] input data for determining individual differences (i.e., a series
            #     of symmetric matrices, number of vars * number of individuals)
            } else if (self$options$mdeMDS == "Ind") {
                private$.prpCfg(crrDim = self$options$dimInd, numRow = length(c(self$options$varInd, self$options$nmeInd)))
            }
        },
        .run = function() {
            crrMDS <- NULL
            # [1] symmetric input data (e.g., distances, correlations, etc.)
            if        (self$options$mdeMDS == "Sym") {
                if ((length(c(self$options$varSym, self$options$nmeSym)) >= self$options$dimSym + 1) &&
                    (nrow(self$data) == length(self$options$varSym))) {
                    crrMDS <- mdsSym(crrDta = self$data,           varSym = self$options$varSym, nmeSym = self$options$nmeSym,
                                     xfmSym = self$options$xfmSym, dimSym = self$options$dimSym, lvlSym = self$options$lvlSym)
                }
            # [2] input data as “usual” data matrix (i.e., vars in columns, units / subjects in rows)
            } else if (self$options$mdeMDS == "Raw") {
                if ((self$options$dirRaw == "col" && length(self$options$varRaw) >= self$options$dimRaw + 1) ||
                    (self$options$dirRaw == "row" && nrow(self$data)             >= self$options$dimRaw + 1)) {
                    crrMDS <- mdsRaw(crrDta = self$data,           varRaw = self$options$varRaw, nmeRaw = self$options$nmeRaw,
                                     xfmRaw = self$options$xfmRaw, dirRaw = self$options$dirRaw, dimRaw = self$options$dimRaw,
                                     lvlRaw = self$options$lvlRaw)
                }
            # [3] input data for determining individual differences (i.e., a series
            #     of symmetric matrices, number of vars * number of individuals)
            } else if (self$options$mdeMDS == "Ind") {
                if ((length(c(self$options$varInd, self$options$nmeInd)) >= self$options$dimInd + 1) &&
                    (nrow(self$data) %% length(self$options$varInd) == 0)) {
                    crrMDS <- mdsInd(crrDta = self$data,           varInd = self$options$varInd, nmeInd = self$options$nmeInd,
                                     xfmInd = self$options$xfmInd, dimInd = self$options$dimInd, lvlInd = self$options$lvlInd)
                }
            }

            if (!is.null(crrMDS)) {
                private$.fllCfg(crrMDS)
                for (nmeFig in c("figCfg", "figHst", "figRes", "figShp", "figStr", "figWgh")) self$results[[nmeFig]]$setState(crrMDS)
            }
            private$.genInf(crrMDS)
        },
        # handling general and additional information messages
        .genInf = function(crrMDS = NULL) {
            crrInf <- self$results$genInf
            crrMde <- self$options$mdeMDS
            if (!is.null(crrMDS)) {
                outInf <- c(sprintf(paste0("<p>Estimated <strong>%s</strong> (of type \"%s\") with %d objects in %d iterations.</p>",
                                           "<p>Stress-1 value: <strong>%.4f</strong></p>"),
                                    crrMDS$model, crrMDS$type, ifelse(is(crrMDS, "smacofR") && self$options$dirRaw == "row", crrMDS$nind, crrMDS$nobj),
                                    crrMDS$niter, crrMDS$stress))
            } else {
                outInf <- c(genMDS, getVar(paste0("gen", crrMde)))
            }
            if (length(outInf) > 0 && all(nchar(outInf) > 0)) {
                crrInf$setContent(paste0("<p>", paste0(outInf, collapse = "</p><p>"), "</p>"))
            }
        },
        .addInf = function() {
             crrInf <- self$results$addInf
             crrMde <- self$options$mdeMDS
             crrInf$setContent(paste0("<p>", paste0(c(addMDS, getVar(paste0("add", crrMde))), collapse = "</p><p>"), "</p>"))
        },
        # prepare and fill configuration table
        .prpCfg = function(crrDim = 2, numRow = 0) {
            if (numRow < 1) return(invisible(NULL))
            crrTbl <- self$results$tblCfg
           
            colNme <- c(sprintf("Dimension %d", seq(crrDim)), rep("SPP", self$options$clmSPP))
            valRow <- stats::setNames(as.list(rep("", length(colNme))), colNme)
            for (i in seq_along(colNme)) crrTbl$addColumn(name = gsub("Dimension ", "D", colNme[i]), title = colNme[i], type = 'number')
            for (i in seq(1, numRow))   crrTbl$addRow(rowKey = i, values = valRow)
        },
        .fllCfg = function(crrMDS = NULL) {
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
            if (is.null(crrMDS)) return(NULL)

            crrCfg <- nmeCfg(crrMDS, self$options$dirRaw)
            crrDta <- cbind(as.data.frame(crrMDS[[crrCfg]]), data.frame(pntSze = bblPnt(crrMDS[[nmeSPP(crrCfg)]], self$options$cfgBbl)))
            crrRng <- c(rndMnM(min(vapply(crrDta[, c("D1", "D2")], min, numeric(1)))),
                        rndMnM(max(vapply(crrDta[, c("D1", "D2")], max, numeric(1)))))
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
            if (is.null(crrMDS)) return(NULL)
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
            if (is.null(crrMDS)) return(NULL)
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
            if (is.null(crrMDS)) return(NULL)
            
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
            if (is.null(crrMDS)) return(NULL)
            
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
            if (is.null(crrMDS)) return(NULL)
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
