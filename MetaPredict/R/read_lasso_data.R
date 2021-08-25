#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' Each file should end in "-ko.tsv" and contain a column of KEGG Orthology terms with column title 'KO'
#' and a column with their associated HMM or Blast E-values titled 'E-value'.
#' @param metadata_file A flatfile with a required column called "filepath". The filepath column must contain the full or relative path to each genome annotation file; Note: files can be from different directories. A second and optional column is genome_name. It should contain the name the user would like to use for each genome. Columns can be in any order, and the rows can also be in any order.
#' @param metadata_df Metadata can optionally be loaded in as a pre-existing dataframe that contains four columns: genome_name, taxonomy, completeness, and filename. Default is FALSE. Can optionally set metadata_df equal to the object name of a pre-existing metadata dataframe.
#' @param kofamscan If the input file or files are output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE.
#' @param cutoff The desired E-value cutoff. Default is 1e-7.
#' @param delim The delimiter of the input files. Default is tab.
#' @importFrom magrittr "%>%"

#' @export
read_data <- function(metadata_file = NULL, metadata_df = FALSE,
                             kofamscan = TRUE, cutoff = 1e-7, delim = '\t', custom = FALSE) {

  cli::cli_h1('Formatting genomic data')

  if (all(!is.null(metadata_file) & metadata_df == FALSE)) {
    metadata_tbl <- readr::read_delim(metadata_file, col_names = TRUE, delim = delim, col_types = readr::cols())

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
    annotations <- purrr::map2(metadata_tbl$filepath, metadata_tbl$genome_name, ~ read_custom(.x, .y, cutoff = cutoff))
  } else if (all(kofamscan == TRUE & custom == TRUE)) {
    cli::cli_alert_danger('Error: Detected multiple TRUE arguments. Please set either the "kofamscan" or "custom" argument equal to TRUE (the other one should be set to FALSE).')
    stop()
  } else {
    cli::cli_alert_danger('Error: something went wrong.')
    stop()
  }

  cli::cli_alert_success('All done.')
  cli::cli_alert_info('Used E-value cutoff: {cutoff}')
  return(annotations)
}


#' @export
read_kofam <- function(.data, .genome_name, cutoff = 1e-7) {
    readr::read_delim(.data, col_types = readr::cols(), delim = '\t') %>%
      dplyr::filter(!stringr::str_detect(`E-value`, '---')) %>%
      dplyr::select(`gene name`, KO, `E-value`) %>%
      dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
      dplyr::filter(`E-value` <= cutoff | `E-value` == 0) %>%
      dplyr::group_by(`gene name`) %>%
      dplyr::filter(`E-value` == min(`E-value`)) %>%
      dplyr::ungroup() %>%
      dplyr::select(KO) %>%
      dplyr::rename(k_number = KO) %>%
      dplyr::mutate(genome_name = .genome_name, .before = 1)
}








# needs to be added to internal data
#filler <- readRDS('training-data-072821/allTrainingFeatureTblKnumberColnames-080421.rda')




read_custom <- function(.data, cutoff = 1e-7, input_type = NULL) {
  .data %>%
    {if (input_type == 'genome_name') anno_type(., cutoff = cutoff, input_type = 'genome_name')
      else if (input_type == 'metagenome_name') anno_type(., cutoff = cutoff, input_type = 'metagenome_name')
      else stop(cli::cli_alert_danger(
        "Error: Names column for input genomes must be 'genome_name', or 'metagenome_name' for input metagenomes column. See usage() for more details"))}
}



