#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' Each file should end in "-ko.tsv" and contain a column of KEGG Orthology terms with column title 'KO'
#' and a column with their associated HMM or Blast E-values titled 'E-value'.
#' @param filePath The full path to the input file or files
#' @param filePattern The file extension pattern. If using multiple files, this will read in all files in the listed directory
#' @param kofamscan If the input file or files are output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @param evalue The desired E-value cutoff. Default is 0.001
#' @param delim The delimiter of the input files. Default is tab
#' @importFrom magrittr "%>%"
#' @export
read_data <- function(filePath, filePattern, kofamscan = TRUE, evalue = 1e-3, delim = '\t') {
  message('Parsing HMM/Blast hits and E-values into MetaPredict. Using E-value cutoff: ', evalue)
  setwd(filePath)
  index <- list.files(path = filePath, pattern = '*-ko.tsv')
  res <- list()

  res = furrr::future_map(1:length(index), .progress = T, ~ {
    col <- sub('(.*)-ko.tsv', '\\1', index[.x], perl = T)
    res[[.x]] <- readr::read_delim(index[.x], col_types = readr::cols(), delim = delim)
    if (kofamscan == TRUE) {
      res[[.x]] <- res[[.x]] %>% dplyr::slice(-1)
    }
    res[[.x]] <- res[[.x]] %>%
      dplyr::select(KO, `E-value`) %>%
      dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
      dplyr::filter(`E-value` <= evalue) %>% #e-value is set by argv, default = 1e-3
      dplyr::select(KO) %>%
      dplyr::rename(!!col := KO) %>%
      dplyr::filter(!(duplicated(.data[[col]])))
    })

  res <- res %>%
    dplyr::tibble() %>%
    tidyr::unnest(cols = dplyr::everything()) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'organism',
                        values_to = 'ko_term', values_drop_na = T) %>%
    dplyr::group_by(organism)

  return(res)
}
