if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("genMDS", "genSym", "genRaw", "genInd", "addSym", "addRaw", "addInd"))
}

# help messages
genMDS <- c(paste("When performing Multidimensional Scaling, you first need to decide what type of data you have.",
                  "jMDS permits three types: (1) symmetric matrices (e.g., correlations, distances), (2) raw data",
                  "(with participants or other units of obsrevation in rows, and variables in columns), and (3)",
                  "matrices that permit to compare individual differences (several symmetric matrices, one per",
                  "person, concatenated over rows)."),
                  "Use the tabs (“Mode”) in the analysis UI to determine the type of input data.")
genSym <-   paste("The symmetric matrices that serve as input for this form of MDS, can have several formats:",
                  "(a) full matrices, (b) triangular matrices, and (c) sparse matrices. If the data are not",
                  "already distances, they need to be converted using the drop-down-menu “Transform Similarities",
                  "to Distances”. See Additional Information at the end of the Results output for more detailed",
                  "explanations.")
genRaw <-   paste("When conducting a MDS with raw data, one first needs to decide what should be compared, rows",
                  "(i.e., persons / units of observation) or columns (i.e., variables). This can be set with the",
                  "radio buttons underneath the variable boxes in the analysis UI. Please note that rows in the",
                  "original data containing NAs are automatically excluded. Typically, the data are not already",
                  "distances. That is, they need to be converted using the drop-down-menu “Transform Similarities to",
                  "Distances”. Depending on which transformation is chosen, the data are either analyzed using",
                  "<code>smacofRect</code> (applies when using the transformations “Reverse” or “Rank”), or using",
                  "<code>smacofSym</code> (applies when using transformations involving correlations or distances).",
                  "See Additional Information at the end of the Results output for more detailed explanations.")
genInd <-   paste("The matrices that serve as input for this form of MDS are matrices (one per individual) that are",
                  "concatenated (over rows). These individual matrices can have one of the following formats: (a)",
                  "full matrices, (b) triangular matrices, and (c) sparse matrices. If the data are not already",
                  "distances, they need to be converted using the drop-down-menu “Transform Similarities to",
                  "Distances”. See Additional Information at the end of the Results output for more detailed",
                  "explanations.")
addSym <- c(      "<h2>Additional Information</h2>",
            paste("<strong>Shape of data matrices</strong>: A <em>full matrix</em> contains both upper and lower",
                  "triangular matrix as well as values in the main diagonal. A <em>triangular matrix</em> contains",
                  "either the upper or the lower triangular matrix, and it can but it doesn't need to contain values",
                  "in the main diagonal. A <em>sparse matrix</em> contains only either the upper or the lower",
                  "triagonal matrix, but has the main diagonal removed. That is, either the first - for upper - or",
                  "the last column - for lower triangulars - are missing. To be able to restore all variable names,",
                  "an additional column with variable names is required."),
            paste("<strong>Type of MDS and ties</strong>: The MDS can either be metric or non-metric (ordinal). For",
                  "“ordinal”, there are different choices how to deal with tied ranks: “primary” - breaks ties, i.e.,",
                  "distances with the same value may be assigned different ranks (with lower ranks assigned to the",
                  "values that appear first in the data); “secondary” - keep ties tied, i.e., distances with the same",
                  "value are assigned the same rank; or “tertiary” - requires that the means of the tie blocks are in",
                  "the correct order."),
            paste("<strong>Tranforming similarities to distances</strong>: If the data are not already distances",
                  "(choose the option “Already Distances (no transform.)” in such case), they need to be tranformed.",
                  "An overview of the mathematical operations behind the different transformations can be found in",
                  "Table 2 of the vignette to the smacof R-package;",
                  "https://cran.r-project.org/web/packages/smacof/vignettes/smacof.pdf)."),
            paste("<strong>Example datasets</strong>: jMDS comes with several example datasets that can be used to",
                  "understand the different possible shapes of data matrices, and to try out the different modes of",
                  "calculating a MDS. These datasets can be accessed using the file menu (☰) - Open - Data Library,",
                  "and then chosing the folder “Multidimensional Scaling (MDS) for jamovi”. Each dataset has tags,",
                  "describing with which mode it should be used (“Symmetric”, “Raw Data”, or “Individual Diff.”), and",
                  "whether the data are already “Distances” or need to be transformed (“Correlations”,",
                  "“Similarities”)."))
addRaw <- c(addSym[1],  # heading is the same as for addSym
            paste("<strong>Shape of data matrices</strong>: Data matrices have the shape that is usual for most / all",
                  "other analyses included in jamovi, namely individuals (or other units of observation) as rows",
                  "and variables as columns."),
            addSym[3],  # description of MDS types and tie handling is the same as for addSym
            paste("<strong>Tranforming similarities to distances</strong>: Typically, raw data are not already",
                  "distances (i.e., ranks within each row, if columns are to be compared). If they are distamces,",
                  "choose “Already Distances (no transform.)”, otherwise choose one of the following tranformations:",
                  "If the data are preferences or ratings, they may be transformed by beeing reversed or ranked.",
                  "When using one of these (or no) tranformation, the shape (rows and columns) of the original data",
                  "is preserved and <code>smacofRect</code> is used for calculating the MDS solution. Alternatively,",
                  "the transformations can be either correlating the data or calculating distances from the raw data.",
                  "If correlations are calculated (over rows or columns, as indicated; and using either parametric -",
                  "Pearson - or non-parametric - Kendall or Spearman - correlations), these correlations are then",
                  "converted into distances. If distance measures are calculated (over rows or columns, as indicated),",
                  "these can be Euclidian, Manhattan, Minkowski or Jaccard distances, and one can further determine",
                  "whether these distance measures shall be calculated on the original or z-standardized data. Using",
                  "correlations or distance measures as transformations changes the shape of the data matrix into a",
                  "symmetric matrix (either of the size rows x rows or columns x columns), which is then analyzed",
                  "using <code>smacofSym</code>."),
            addSym[5]) # comment about the example datasets is the same as for addSym
addInd <- c(addSym[1],  # heading is the same as for addSym
            paste(addSym[2], "When comparing individual differences, each individual is described by one matrix of",
                  "these matrix types. The matrices are concatenated over rows. That is, if V is the number of",
                  "variables assigned to “Variables for MDS”, the first V rows are the data for the first individual,",
                  "the second V rows for the second individual, and so on."),
            addSym[3], # description of MDS types and tie handling is the same as for addSym
            addSym[4], # description of transformation operations is the same as for addSym
            addSym[5]) # comment about the example datasets is the same as for addSym
