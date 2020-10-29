#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' Each file should end in "-ko.tsv" and contain a column of KEGG Orthology terms with column title 'KO'
#' and a column with their associated HMM or Blast E-values titled 'E-value'.
#' @param Data A dataframe or list of dataframes of HMM/Blast results which include KEGG Orthology terms and E-values in two of the columns of each dataframe
#' @param File The full or relative path to a single input flatfile
#' @param filePath The full path to the input file or files
#' @param filePattern The file extension pattern. If using multiple files, this will read in all files in the listed directory
#' @param kofamscan If the input file or files are output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @param evalue The desired E-value cutoff. Default is 0.001
#' @param delim The delimiter of the input files. Default is tab
#' @importFrom magrittr "%>%"
#' @export
new_read_genome_data <- function(Data = NULL, File = NULL, filePath = NULL, filePattern = '-ko.tsv', #still needs to accept single input files more easily...
                             kofamscan = TRUE, evalue = 1e-3, delim = '\t') {

  message('Parsing HMM/Blast hits and E-values into MetaPredict. Using E-value cutoff: ', evalue)
  res <- list()

  if (is.null(Data) & is.null(File)) {
    files <- list.files(path = filePath, pattern = filePattern, full.names = TRUE)
    index <- length(files)
  }
  if (is.null(Data) & !is.null(File)) {
    index <- length(File)
  }
  if (is.list(Data) & purrr::every(Data, ~ {
    grepl(paste(names(.x), collapse = '|'),
          c('KO', 'evalue|e-value'), ignore.case = T)
  })) {
    index <- length(Data)
    orgName <- substitute(Data)
  }
  res <- purrr::map(1:index, .progress = T, ~ {
    if (is.null(Data) & is.null(File)) {
      col <- sub(paste('.*\\/(.*)', filePattern, sep = ''), '\\1', files[.x], perl = T)
      res[[.x]] <- readr::read_delim(files[.x], col_types = readr::cols(), delim = delim)

    } else if (is.null(Data) & !is.null(File)) {
      col <- sub('.*\\/(.*)\\..*', '\\1', File, perl = T)
      res[[.x]] <- readr::read_delim(File, col_types = readr::cols(), delim = delim)

    } else if (is.data.frame(Data) & length(list(Data)) == 1) {
      res[[.x]] <- Data

    } else if (is.list(Data) & length(Data) >= 1) {
      res[[.x]] <- Data[[.x]]

    } else {
      stop('Data must be of class dataframe, list of dataframes, or an input flatfile/flatfile(s).')
    }

    if (kofamscan == TRUE) {
      res[[.x]] <- res[[.x]] %>% dplyr::slice(-1)
    }
    res[[.x]] <- res[[.x]] %>%
      dplyr::select(`gene name`, KO, `E-value`) %>% #add in gene_name info, and keep E-values, don't filter them!
      dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
      dplyr::filter(`E-value` <= evalue) %>%
      dplyr::select(KO, `gene name`) %>%
      dplyr::rename(Gene = `gene name`)

    if (is.null(Data)) {
      res[[.x]] <- res[[.x]] %>%
        dplyr::filter(!(duplicated(KO))) %>%
        dplyr::summarize(ko_term = paste0(KO, collapse = ' '),
                         Gene = paste0(Gene, collapse = ' ')) %>%
        dplyr::mutate(organism := !!col)

    } else if (is.data.frame(Data) & length(list(Data)) == 1) {
      res[[.x]] <- res[[.x]] %>%
        dplyr::filter(!(duplicated(KO))) %>%
        dplyr::summarize(ko_term = paste0(KO, collapse = ' '),
                         Gene = paste0(Gene, collapse = ' ')) %>%
        dplyr::mutate(organism := !!deparse(orgName))

    } else {
      res[[.x]] <- res[[.x]] %>%
        dplyr::filter(!(duplicated(KO))) %>%
        dplyr::summarize(ko_term = paste0(KO, collapse = ' '),
                         Gene = paste0(Gene, collapse = ' ')) %>%
        dplyr::mutate(organism := !!deparse(paste(orgName, .x, sep = '_')))
    }
  })
  res <- res %>%
    dplyr::tibble() %>%
    tidyr::unnest(cols = dplyr::everything()) %>%
    dplyr::group_by(organism)

  return(res)
}
