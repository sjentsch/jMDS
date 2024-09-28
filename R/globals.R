if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("genMDS", "genSym", "genRaw", "genInd", "addMDS", "addSym", "addRaw", "addInd"))
}

# help messages
genMDS <- c(paste("When performing Multidimensional Scaling, you first need to decide what type of data you have.",
                  "jMDS permits three types: (1) symmetric matrices (e.g., correlations, distances), (2) raw data",
                  "(with participants or other units of obsrevation in rows, and variables in columns), and (3)",
                  "matrices that permit to compare individual differences (several symmetric matrices, one per",
                  "person, concatenated over rows)."),
                  "Use the tabs (\"Mode\") in the analysis UI to determine the type of input data.")
genSym <- c(paste("The symmetric matrices that serve as input for this form of MDS, can have several formats:",
                  "(a) full matrices (containing both upper and lower triangular matrix and the main diagonal),",
                  "(b) triangular matrices (containing either the upper or the lower triangular matrix as well",
                  "as the main diagonal), and (c) sparse matrices (containing only either the upper or the lower",
                  "triangular matrix; NB: when using sparse matrices, they have to contain distances)."),
            paste("If the data in the symmetric matrix are not already distances, they need to be converted.",
                  "Possible transformations can be found in the drop-down-menu \"Transform Similarities to",
                  "Distances\" (for details what these transformations do, see Table 2 in the vignette to the",
                  "smacof R-package; https://cran.r-project.org/web/packages/smacof/vignettes/smacof.pdf)."))
genRaw <- c(paste("When conducting a MDS with raw data, one first needs to decide what should be compared, rows",
                  "(i.e., persons / units of observation) or columns (i.e., variables). This can be set with the",
                  "radio buttons underneath the variable boxes in the analysis UI."),
            paste("The next decision is (a) whether the data are already distances (i.e., ranks within each row",
                  "if columns are to be compared), (b) whether the data are preferences or ratings that need to",
                  "be reversed, (c) whether correlations are calculated (over rows or columns, as indicated; ",
                  "either parametric - Pearson - or non-parametric - Kendall or Spearman) and then converted into",
                  "distances, or (d) whether distance measures should be calculated (over rows or columns, as",
                  "indicated; Euclidian, Manhattan, Minkowski or Jaccard) on the original or z-standardized data.",
                  "(a) and (b) preserve the shape (rows and columns) of the original data and smacofRect is used",
                  "for the analysis, (c) and (d) result in symmetric matrices that are analyzed using smacofSym."),
                  "Please note that rows in the original data containing NAs are automatically excluded.")
genInd <- c(paste("Gen Ind 1.1",
                  "Gen Ind 1.2"),
                  "Gen Ind 2")
addMDS <- c(paste("Add MDS 1.1",
                  "Add MDS 1.2"),
                  "Add MDS 2")
addSym <- c(paste("Details Sym 1.1",
                  "Details Sym 1.2"),
            paste("Details Sym 2.1",
                  "Details Sym 2.2"))
addRaw <- c(paste("Details Raw 1.1",
                  "Details Raw 1.2"),
            paste("Details Raw 2.1",
                  "Details Raw 2.2"))
addInd <- c(paste("Details Ind 1.1",
                  "Details Ind 1.2"),
            paste("Details Ind 2.1",
                  "Details Ind 2.2"))
