#' Title
#'
#' @param taxa A table containing taxa, either from a `download`, a `download_list` or
#' @param output A string representing the file path for export.  Can be NULL.
#'
#' @return A \code{data.frame}.
#' @export
#'
#' @examples
#'

generate_tables <- function(taxa, output = NULL) {

  if (!any(c('download', 'download_list') %in% class(taxa))) {
    assertthat::assert_that('data.frame' %in% class(taxa), msg = "taxa must be a neotoma download object or a data.frame.")
    taxon <- data.frame(target = colnames(taxa), match = NA, stringsAsFactors = FALSE)

  } else {

    taxon <- data.frame(neotoma::taxa(taxa),
                        match = NA,
                        stringsAsFactors = FALSE)
  }

  if (is.null(output)) {
    return(taxon)
  } else {
    write.csv(taxon, output)
  }
}
