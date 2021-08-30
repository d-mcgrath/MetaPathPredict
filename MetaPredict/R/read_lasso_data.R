#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' @param input_dir Character. The path to the directory containing the input annotation files to use, if a metadata file or dataframe is not supplied. Default is NULL.
#' @param pattern Character. If an input_dir argument is supplied, this is the regex pattern to use to locate the files within the directory containing the annotation files. Default is '*.tsv'.
#' @param delim Character. The delimiter of the input files. Default is tab.
#' @param metadata_file Character. The path to a flatfile with a required column called "filepath". The filepath column must contain the full or relative path to each genome annotation file; Note: files can be from different directories. A second and optional column is genome_name. It should contain the name the user would like to use for each genome. Columns can be in any order, and the rows can also be in any order. Default genome names will be provided if the genome_name column is not present in metadata_file.
#' @param metadata_df Dataframe or Tibble. Metadata can optionally be provided as a pre-existing dataframe that contains the column "filepath". A second and optional column is genome_name. It should contain the name the user would like to use for each genome. Columns can be in any order, and the rows can also be in any order. Default genome names will be provided if the genome_name column is not present in metadata_df.
#' @param kofamscan Logical. If the input file or files are output from the Kofamscan command line tool, leave this argument set to TRUE, otherwise set it to FALSE and set the custom argument to TRUE. Default is TRUE.
#' @param cutoff Numeric. The desired E-value or score cutoff. Default is 1e-7.
#' @param custom Logical. If input annotations are not the output of Kofamscan, the column containing KEGG K numbers must be called "k_number". If annotations contain scores - HMM E-values, blast bitscores, etc., then the gene name of each gene must be in a column "gene_name", and the column with the score must be called "evalue" for E-values from HMM hits, or "score" for any other score, e.g., a bitscore.
#' @param score_type Character: "evalue", "score", or "none". If "evalue", there must be a column in input annotation files called "evalue" for E-values from HMM hits, or "score" for any other score, e.g., a bitscore. If "none", no score column for annotations is required.
#' @importFrom magrittr "%>%"

#' @export
read_data <- function(input_dir = NULL, pattern = '*.tsv', delim = '\t', metadata_file = NULL, metadata_df = NULL,
                      kofamscan = TRUE, cutoff = 1e-7, custom = FALSE, score_type = c('evalue', 'score', 'none')) {

  cli::cli_h1('Formatting data')

  score_type <- match.arg(score_type)

  if (!is.null(input_dir)) {
    if (dir.exists(input_dir)) {
      files <- list.files(path = input_dir, pattern = pattern, full.names = TRUE)
      genome_names <- list.files(path = input_dir, pattern = pattern, full.names = FALSE)
      annotations <- read_from_dir(files, genome_names, cutoff = cutoff, delim = delim,
                                   kofamscan = kofamscan, custom = custom, score_type = score_type)
    } else {
      cli::cli_alert_danger('Error: Directory given for "input_dir" argument does not exist.')
    }
  } else {

    if (all(!is.null(metadata_file) & is.null(metadata_df))) {
      metadata_tbl <- readr::read_delim(metadata_file, delim = delim, col_types = readr::cols())

    } else if (all(is.null(metadata_file) & is.data.frame(metadata_df))) {
      metadata_tbl <- metadata_df

    } else if (all(!is.null(metadata_file) & is.data.frame(metadata_df))) {
      cli::cli_alert_danger('Error: either a metadata filepath should be provided to the metadata_file argument, or a dataframe object should be provided to the metadata_df argument. Inputs were detected for both the metadata_file and metadata_df arguments; please select only one and remove the other one.')
      stop()
    } else {
      cli::cli_alert_danger('Error: no input was detected to either of the metadata_file/metadata_df arguments. Please make sure to supply metadata input to one of these arguments.')
      stop()
    }

    if ('filepath' %in% colnames(metadata_tbl)) {
      if ('genome_name' %in% colnames(metadata_tbl)) {
        if (any(duplicated(metadata_tbl$genome_name))) {
          cli::cli_alert_danger('Error: please make sure that all names in the genome_name column are unique.')
          stop()
        }
        if (any(is.na(metadata_tbl$genome_name))) {
          cli::cli_alert_warning('Warning: detected ', length(which(is.na(metadata_tbl$genome_name))), 'NA values in the genome_name column. Adding in default genome names for these missing values.')
          metadata_tbl <- metadata_tbl %>%
            dplyr::mutate(genome_name = dplyr::case_when(is.na(genome_name) ~ paste0('genome_', dplyr::cur_group_rows()),
                                                         TRUE ~ genome_name))
        }
      } else {
        cli::cli_alert_info('genome_name column not detected. Adding genome_name column to metadata, using default name values.')
        metadata_tbl <- metadata_tbl %>%
          dplyr::mutate(genome_name = paste0('genome_', seq_along(filepath)))
      }
    } else {
      potential_cases <- metadata_tbl %>%
        dplyr::select_if(~ any(stringr::str_detect(.x, stringr::regex('\\/|\\~|\\.')))) %>%
        colnames()
      cli::cli_alert_danger('Error: Does your metadata file/dataframe contain the column name: "filepath" with the full filepath to each of the genome annotations? If not, please add a column named "filepath" to your metadata file/dataframe that contains the genome annotation full filepaths.')
      if (potential_cases > 0) {
        message('Do any of the following columns contain your genome annotation filepaths: ', potential_cases, collapse = '; ', '\n', 'If yes, rename the column containing the filepaths to "filepath". Otherwise, please add a column named "filepath" to your metadata file that contains the genome annotation full filepaths.')
      }
      stop()
    }

    if (all(kofamscan == TRUE & custom == FALSE)) {
      annotations <- purrr::map2(metadata_tbl$filepath, metadata_tbl$genome_name, ~ read_kofam(.x, .y, cutoff = cutoff))
    } else if (all(kofamscan == FALSE & custom == TRUE)) {
      annotations <- purrr::map2(metadata_tbl$filepath, metadata_tbl$genome_name, ~ read_custom(.x, .y, cutoff = cutoff, score_type = score_type))
    } else if (all(kofamscan == TRUE & custom == TRUE)) {
      cli::cli_alert_danger('Error: Detected multiple TRUE arguments. Please set either the "kofamscan" or "custom" argument equal to TRUE (the other one should be set to FALSE).')
      stop()
    } else {
      cli::cli_alert_danger('Error: something went wrong.')
      stop()
    }
  }
  cli::cli_alert_success('All done.')
  if (score_type == 'evalue' | score_type == 'score') {
    cli::cli_alert_info('Used annotation scoring cutoff: {cutoff}')
  } else {
    cli::cli_alert_info('Did not filter input annotations by any scoring metric.')
  }
  return(annotations)
}



#' @export
read_from_dir <- function(.files, .genome_names, cutoff = cutoff, delim = delim,
                          kofamscan = kofamscan, custom = custom, score_type = score_type) {
  if (all(kofamscan == TRUE & custom == FALSE)) {
    annotations <- purrr::map2(.files, .genome_names, ~ read_kofam(.x, .y, cutoff = cutoff))
  } else if (all(kofamscan == FALSE & custom == TRUE)) {
    annotations <- purrr::map2(.files, .genome_names, ~ read_custom(.x, .y, cutoff = cutoff, delim = delim, score_type = score_type))
  } else if (all(kofamscan == TRUE & custom == TRUE)) {
    cli::cli_alert_danger('Error: Detected multiple TRUE arguments. Please set either the "kofamscan" or "custom" argument equal to TRUE (the other one should be set to FALSE).')
    stop()
  } else {
    cli::cli_alert_danger('Error: something went wrong.')
    stop()
  }
  return(annotations)
}



#' @export
read_kofam <- function(.data, .genome_name, cutoff = 1e-7) {
  readr::read_tsv(.data, col_types = readr::cols()) %>%
    dtplyr::lazy_dt() %>%
    dplyr::filter(!stringr::str_detect(`E-value`, '---')) %>%
    dplyr::select(`gene name`, KO, `E-value`) %>%
    dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
    dplyr::filter(`E-value` <= cutoff | `E-value` == 0) %>%
    dplyr::group_by(`gene name`) %>%
    dplyr::filter(`E-value` == min(`E-value`)) %>%
    dplyr::ungroup() %>%
    dplyr::select(KO) %>%
    dplyr::rename(k_number = KO) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(genome_name = .genome_name, .before = 1)
}





## needs to be updated
read_custom <- function(.data, .genome_name, cutoff = 1e-7, delim = delim, score_type = 'evalue') {

  if (score_type == 'evalue') {
    result <- readr::read_delim(.data, col_types = readr::cols(), delim = delim) %>%
      dtplyr::lazy_dt() %>%
      dplyr::select(gene_name, k_number, dplyr::matches('^e_value$|^evalue$|^e-value$')) %>%
      dplyr::rename(e_value = 3) %>%
      dplyr::mutate(e_value = as.numeric(e_value)) %>%
      dplyr::filter(e_value <= cutoff | e_value == 0) %>%
      dplyr::group_by(gene_name) %>%
      dplyr::filter(e_value == min(e_value)) %>%
      dplyr::ungroup() %>%
      dplyr::select(k_number) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(genome_name = .genome_name, .before = 1)
  } else if (score_type == 'score') {
    result <- readr::read_delim(.data, col_types = readr::cols(), delim = delim) %>%
      dtplyr::lazy_dt() %>%
      dplyr::select(gene_name, k_number, score) %>%
      dplyr::mutate(score = as.numeric(score)) %>%
      dplyr::filter(score >= cutoff) %>%
      dplyr::group_by(gene_name) %>%
      dplyr::filter(score == max(score)) %>%
      dplyr::ungroup() %>%
      dplyr::select(k_number) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(genome_name = .genome_name, .before = 1)
  } else if (score_type == 'none') {
    result <- readr::read_delim(.data, col_types = readr::cols(), delim = delim) %>%
      dtplyr::lazy_dt() %>%
      dplyr::select(k_number) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(genome_name = .genome_name, .before = 1)
  } else {
    cli::cli_alert_danger('Error: Issue with score argument Please make sure it is "evalue", "score", or "none". If score argument is "evalue" or "score", please make sure the score column is named "evalue" or "score".')
  }
    return(result)
}
