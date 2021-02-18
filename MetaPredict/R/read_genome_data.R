#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' Each file should end in "-ko.tsv" and contain a column of KEGG Orthology terms with column title 'KO'
#' and a column with their associated HMM or Blast E-values titled 'E-value'.
#' @param Data A dataframe or list of dataframes of HMM/Blast results which include KEGG Orthology terms and E-values in two of the columns of each dataframe, or the full or relative path to a single or multiple input flatfiles containing this information
#' @param taxonCols A dataframe with one column that lists the taxonomy of each input genome, and a second column that lists the estimated genome completeness of each input genome. Optionally, can be three columns where column 1 is the genome name, column two is the genome taxonomy, and column 3 is the estimated genome completeness of each input genome. Order of the columns should match the order of the input genomes.
#' @param filePath The full path to the input file or files.
#' @param filePattern The file extension pattern. If using multiple files, this will read in all files in the listed directory.
#' @param kofamscan If the input file or files are output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE.
#' @param evalue The desired E-value cutoff. Default is 0.001.
#' @param delim The delimiter of the input files. Default is tab.
#' @importFrom magrittr "%>%"

#' @export
read_genome_data <- function(Data = NULL, taxonCols = NULL, filePath = NULL, filePattern = '-ko.tsv',
                             kofamscan = TRUE, dram = FALSE, cutoff = 1e-3, delim = '\t') {

  cli::cli_h1('Reading genomic data into MetaPredict')

  if (!is.null(taxonCols)) {
    suppressWarnings(taxon_tibble <- readr::read_delim(taxonCols,
                                                       col_names = c('taxonomy', 'est_comp'), delim = delim,
                                                       col_types = readr::cols())) # NEEDS UPDATE; col length can vary between 2 and 3...
  } else {
    stop(cli::cli_alert_danger('Taxon and estimated genome completness information required. Please provide a flatfile consisting of a column of taxon names of each input genome, and a column of the estimated completeness of each genome, in the same order as the genomes input into MetaPredict.'))
  }
  index <- get_index(.data = Data, filePath = filePath, filePattern = filePattern)

  if (is.null(Data)) {
    files <- list.files(path = filePath, pattern = filePattern, full.names = TRUE)
  }
  if (is.list(Data)) { #colname checks for all forms of Data: NEEDS UPDATE to include all columns that will be REQUIRED
    if (!purrr::every(Data, ~ {
      grepl(paste(colnames(.x), collapse = '|'), c('KO', 'evalue|e-value'), ignore.case = T)})) {
      cli::cli_alert_danger('Error: Colnames  KO and e-value were not detected. make sure input columns are properly labelled.')
      stop()
    } else {.}
  }

  formatted_genomes <- list()
  formatted_genomes <- purrr::map(1:index, .progress = T, ~ { #for each element in the index, index is just number vector- # of inputs
    if (is.list(Data) & length(Data) >= 1) { # MAKE SURE THIS ALWAYS RETURNS 1 FOR DATAFRAMES
      Data <- Data[[.x]]
      orgName = names(Data[.x]) # this is supposed to give a default name each genome tibble in a list if none are provided
    }
    if (!('genome_name' %in% colnames(taxon_tibble))) {
      if (is.null(Data)) {
        files_x <- files[.x]
        orgName <- get_fileName.g(filePattern = filePattern, files = files_x)
      } else if (is.character(Data)) {
        orgName <- get_fileName.g(.data = Data)
      } else if (is.list(Data) & length(Data) >= 1 & !is.data.frame(Data)) {
        orgName <- names(Data[.x])
      } else if (is.list(Data) & is.data.frame(Data) & length(Data) >= 1) {
        orgName <- substitute(Data)
      } else {
        cli::cli_alert_danger('Error: Data must be a dataframe, list of dataframes, a flatfile, or multiple flatfiles. Please see usage() for more information.')
        stop()
      }
    } else {
      orgName <- NULL
    }

    if (is.null(Data)) {
      files_x <- files[.x]
      formatted_genomes[[.x]] <- import_data.g(.data = Data, files = files_x, delim = delim)
    } else {
      formatted_genomes[[.x]] <- import_data.g(.data = Data, delim = delim)
    }

    formatted_genomes[[.x]] <- formatted_genomes[[.x]] %>%
      {if (kofamscan == TRUE) tidy_kofam(., cutoff = cutoff, `E-value`, KO, `gene name`)
        else if (dram == TRUE) tidy_dram(.)
        else if (kofamscan == FALSE & dram == FALSE) tidy_custom_anno(., cutoff = cutoff, input_type = 'genome_name')
        else stop(cli::cli_alert_danger('Error: Issue importing user data. Please see usage() for data import guidelines.'))} %>%
      concatenate(orgName = orgName, taxon = taxon_tibble$taxonomy[[.x]]) %>%
      dplyr::mutate(p_j = taxon_tibble$est_comp[.x])
  }) %>%
    dplyr::tibble() %>%
    tidyr::unnest(cols = dplyr::everything()) %>%
    dplyr::group_by(genome_name) %>%
    dplyr::mutate(data_type = 'genome')

  cli::cli_alert_success('Parsed HMM/Blast hits and E-values')
  cli::cli_alert_info('Used E-value cutoff: {cutoff}')

  return(formatted_genomes)
}



get_index <- function(.data, filePath, filePattern) {
  if (is.character(.data)) {
    index <- length(.data)
  } else if (is.null(.data)) {
    index <- length(list.files(path = filePath, pattern = filePattern))
  } else if (is.list(.data) & purrr::every(.data, ~ {grepl(paste(names(.x), collapse = '|'),
                                                           c('KO', `E-value`), ignore.case = TRUE)})) {
    index <- length(.data)
  } else {cli::cli_alert_danger("Error: Issue generating index for data import.")
    stop()
  }
  return(index)
}



import_data.g <- function(.data, files = NULL, delim) {
  if (is.null(.data)) {
    import <- readr::read_delim(files, col_types = readr::cols(), delim = delim)
  } else if (is.character(.data)) {
    import <- readr::read_delim(.data, col_types = readr::cols(), delim = delim)
  } else if (is.data.frame(Data) & length(list(.data)) == 1) {
    import <- .data
  } else if (is.list(.data) & length(.data) >= 1) {
    import <- .data
  } else {cli::cli_alert_danger(
    'Data input failed. Data must be a dataframe, list of dataframes, or an input flatfile/flatfile(s).')
    stop()
  }
  return(import)
}



get_fileName.g <- function(.data = NULL, filePattern = NULL, files = NULL) {
  if (is.null(.data)) {fileName <- sub(paste('.*\\/(.*)', filePattern, sep = ''), '\\1', files, perl = T) # files was 'files[.x]'
  } else if (is.character(.data)) {fileName <- sub('.*\\/(.*)\\..*', '\\1', .data, perl = T) #is.character() may not work as expected here
  #} else if (is.list(.data)) {fileName <- names(.data)
  } else {cli::cli_alert_danger('Data input failed. Data must be a dataframe, list of dataframes, or an input flatfile/flatfile(s).')
    stop()}
  return(fileName)
}



tidy_kofam <- function(.data, cutoff = 1e-3, ...) {
  .data %>%
    filter(!str_detect(`E-value`, '-----')) %>%
    {if (all(c('E-value', 'KO', 'gene name') %in% colnames(.))) dplyr::select_(., .dots = lazyeval::lazy_dots(...))
      else stop(cli::cli_alert_danger(
        "Error: Columns 'E-value', 'KO', and 'gene name' not detected. These columns are required to read in Kofamscan output files."))} %>%
    dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
    dplyr::filter(`E-value` <= cutoff | `E-value` == 0) %>%
    dplyr::select(KO, `gene name`) %>%
    dplyr::rename(gene_name = `gene name`, k_number = KO)
}



tidy_dram <- function(.data) {
  .data %>%
    {if (all(c('fasta', 'scaffold', 'gene_position',
               'kegg_id', 'kegg_hit') %in% colnames(.)))  dplyr::select(fasta, scaffold, gene_position, kegg_id, kegg_hit)
      else if (all(c('fasta', 'kegg_id', 'kegg_hit') %in% colnames(.))) dplyr::select(fasta, kegg_id, kegg_hit)
      else if (all(c('kegg_id', 'kegg_hit') %in% colnames(.))) dplyr::select(kegg_id, kegg_hit)
      else stop(cli::cli_alert_danger(
        "Error: Columns 'kegg_id' and 'kegg_hit' not detected. These columns are required to read in Dram output files."))} %>%
    dplyr::filter(!is.na(kegg_id), !duplicated(kegg_id))
}



anno_type <- function(.data, cutoff = 1e-3, input_type = NULL) {
  .data %>%
    {if (all(c('k_number', 'gene_name', input_type, 'hit_score') %in% colnames(.)) |
         all(c('k_number', 'gene_name', input_type) %in% colnames(.)) |
         all(c('k_number', 'gene_name') %in% colnames(.)) |
         'k_number' %in% colnames(.))
      dplyr::select(., if (all(c('k_number', 'gene_name', input_type, 'hit_score') %in% colnames(.))) {c(genome_name, gene_name, k_number, hit_score)}
                    else if (all(c('k_number', 'gene_name', input_type) %in% colnames(.))) {c(genome_name, gene_name, k_number)}
                    else if (all(c('k_number', 'gene_name') %in% colnames(.))) {c(k_number, gene_name)}
                    else if ('k_number' %in% colnames(.)) {k_number}
                    else stop(cli::cli_alert_danger("Error: Column 'k_number' was not detected. k_number column is required. [gene_name, name, hit_score columns are optional]"))) %>%
        {if ('hit_score' %in% colnames(.)) dplyr::mutate(., hit_score = as.numeric(hit_score)) %>%
            dplyr::filter(., hit_score <= cutoff | hit_score == 0)
          else (.)}
      else stop(cli::cli_alert_danger('Something went wrong importing custom data. Please see usage().'))
    }
}



tidy_custom_anno <- function(.data, cutoff = 1e-3, input_type = NULL) {
  .data %>%
    {if (input_type == 'genome_name') anno_type(., cutoff = cutoff, input_type = 'genome_name')
      else if (input_type == 'metagenome_name') anno_type(., cutoff = cutoff, input_type = 'metagenome_name')
      else stop(cli::cli_alert_danger(
        "Error: Names column for input genomes must be 'genome_name', or 'metagenome_name' for input metagenomes column. See usage() for more details"))}
}



concatenate <- function(.data, orgName = NULL, taxon = NA, default_name = NULL) {
  .data %>%
    {if (is.null(.data) | is.character(.data)) dplyr::filter(., !(duplicated(k_number))) %>%
        dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
                         gene_name = paste0(gene_name, collapse = ' ')) %>%
        {if (!is.null(orgName)) dplyr::mutate(., genome_name = orgName) #genome_name := !!orgName)
          else (.)} %>%
        dplyr::mutate(., taxonomy = taxon)

      else if (is.data.frame(.data) & length(list(.data)) == 1 & length(.data) >= 1) {
        dplyr::filter(., !(duplicated(k_number))) %>%
          dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
                           gene_name = paste0(gene_name, collapse = ' ')) %>%
          {if (!is.null(orgName)) dplyr::mutate(., genome_name = orgName)  #:= !!deparse(orgName))
            else (.)} %>%
          dplyr::mutate(taxonomy = taxon)

      } else {dplyr::filter(., !(duplicated(k_number))) %>% #THIS NEEDS REWORK
          dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
                           gene_name = paste0(gene_name, collapse = ' ')) %>%
          {if (!is.null(orgName)) dplyr::mutate(., genome_name := !!deparse(paste(orgName, .x, sep = '_'))) # this will fail. needs rework.
            else (.)} %>%
          dplyr::mutate(taxonomy = taxon)}
    }
}
