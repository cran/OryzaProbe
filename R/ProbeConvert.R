#' @title probeConvert
#' @description Convert rice microarray probe ID to RAP-DB ID and merge multiple probes for one gene
#' @param exprMatrix A expression matrix is required. Row names of the matrix should be the microarray probe ID.
#' @param probeMerge A instruction is needed. If probeMerge = F, the multiple probes for one gene will not be merged. If probeMerge = T, the multiple probes for one gene will be merged. The default value is F.
#' @param mergeBy A method for merging multiple probes for one gene if probeMerge = T. Mean, max and min value is available.
#' @return A matrix.
#' @importFrom stats aggregate
#' @examples
#' \donttest{
#' expr <- system.file("test_file", header = TRUE, sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' formatted_expr <- probeConvert(exprMatrix = expr)
#' formatted_expr1 <- probeConvert(exprMatrix = expr, probeMerge = TRUE)
#' }
#' @export


probeConvert <- function(exprMatrix, probeMerge = FALSE, mergeBy = 'mean'){
  if (12 %in% rownames(exprMatrix)){
    OryzaGPL = "GPL6864"
    dataID = GPL6864
    Info = list(OryzaGPL, dataID)
    message("GPL platform is GPL6864")
    mode = 1}
  else if ('Os.1.1.S1_s_at' %in% rownames(exprMatrix)){
    OryzaGPL = "GPL2025"
    dataID = GPL2025
    Info = list(OryzaGPL, dataID)
    message("GPL platform is GPL2025")
    mode = 1}
  else if ('Os01g0100100|COMBINER_EST|CI448596|0' %in% rownames(exprMatrix)){
    OryzaGPL = "GPL8852"
    dataID = GPL8852
    Info = list(OryzaGPL, dataID)
    message("GPL platform is GPL8852")
    mode = 1}
  else {
    message("The platform may not be supported now. GPL2025, GPL6864, and GPL8852 are available. Please check your GPL platform.")
    mode = 0}
  if (mode == 1){
    if (probeMerge == FALSE && Info[1] == 'GPL8852'){
      rowNames <- rownames(exprMatrix)
      ACC.GPL8852 <- sub("\\|.+", "", rowNames)
      ACC.GPL8852.uniq <- make.unique(ACC.GPL8852, sep = ".")
      rownames(exprMatrix) <- ACC.GPL8852.uniq
      return(exprMatrix)}
    else if (probeMerge == FALSE && Info[1] != 'GPL8852'){
      ACC <- Info[[2]][["ACC"]]
      annoID <- Info[[2]][["ID"]]
      exprID <- rownames(exprMatrix)
      num <- which(annoID %in% exprID)
      final_ACC <- ACC[num]
      final_ACC.uniq <- make.unique(final_ACC, sep = ".")
      rownames(exprMatrix) <- final_ACC.uniq
      return(exprMatrix)}
    else if (probeMerge == TRUE && Info[1] == 'GPL8852'){
      rowNames <- rownames(exprMatrix)
      ACC.GPL8852 <- sub("\\|.+", "", rowNames)
      exprMatrix <- as.data.frame(exprMatrix, as.numeric)
      exprMatrix <- cbind(exprMatrix, ACC.GPL8852)
      exprMatrix <- aggregate(x = exprMatrix[,1:ncol(exprMatrix) - 1], by = list(ACC.GPL8852), FUN = mergeBy)
      rownames(exprMatrix) = exprMatrix[, 1]
      exprMatrix <- exprMatrix[, -1]
      return(exprMatrix)}
    else if (probeMerge == TRUE && Info[1] != 'GPL8852'){
      ACC <- Info[[2]][["ACC"]]
      annoID <- Info[[2]][["ID"]]
      exprID <- rownames(exprMatrix)
      exprMatrix <- as.data.frame(exprMatrix, as.numeric)
      exprMatrix <- cbind(exprMatrix, exprID)
      exprMatrix <- aggregate(x = exprMatrix[,1:ncol(exprMatrix) - 1], by = list(exprID), FUN = mergeBy)
      rownames(exprMatrix) = exprMatrix[, 1]
      exprMatrix <- exprMatrix[, -1]
      return(exprMatrix)}}
}
