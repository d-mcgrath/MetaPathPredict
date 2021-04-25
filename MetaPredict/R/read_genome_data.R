#' This function is used to read in user data in the form of a one or more flatfiles, default is TSV format.
#' Each file should end in "-ko.tsv" and contain a column of KEGG Orthology terms with column title 'KO'
#' and a column with their associated HMM or Blast E-values titled 'E-value'.
#' @param Data A dataframe or list of dataframes of HMM/Blast results which include KEGG Orthology terms and E-values in two of the columns of each dataframe, or the full or relative path to a single or multiple input flatfiles containing this information
#' @param metadata A flatfile with one column "taxonomy" that lists the taxonomy of each input genome, and a second column "completeness" that lists the estimated genome completeness of each input genome. The completeness column should be in the format of percent completeness - e.g., a number between 0 and 100 for each genome. Taxonomy should be the NCBI taxonomic name for each genome, from Domain to Species level - it can optionally be "Unknown" if the taxonomic name is not known. The genome_name column should contain the name the user would like to use for each genome. The filename column must contain the full or relative path to each genome annotation file; Note: files can be from different directories.  If metaColnames is TRUE, metadata should have four columns including genome_name, taxonomy, completeness, and filename for each input genome; information for each genome should occupy a single row of the flatfile. Columns can be in any order, and the rows can also be in any order. NOTE: if metaColnames is FALSE - which is not recommended - three unnamed columns will be required that contain the taxonomy, completeness, and filename for each genome.
#' @param metaColnames Logcial. If TRUE, the metadata flatfile contains the following columns in any order: genome_name, taxonomy, completeness, filename. If FALSE, the flatfile contains THREE unnamed columns in any order that contain the taxonomy, estimated completeness, and filename for each genome; it is not recommended to use this setting if possible.
#' @param metadata_df Metadata can optionally be loaded in as a pre-existing dataframe that contains four columns: genome_name, taxonomy, completeness, and filename. Default is FALSE. Can optionally set metadata_df equal to the object name of a pre-existing metadata dataframe.
#' @param kofamscan If the input file or files are output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE.
#' @param evalue The desired E-value cutoff. Default is 0.001.
#' @param delim The delimiter of the input files. Default is tab.
#' @importFrom magrittr "%>%"

#' @export
read_genome_data <- function(Data = NULL, metadata = NULL, metaColnames = TRUE, metadata_df = FALSE, #filePath = NULL, filePattern = '-ko.tsv',
                             kofamscan = TRUE, dram = FALSE, cutoff = 1e-5, delim = '\t') {

  cli::cli_h1('Formatting genomic data')

  if (all(!is.null(metadata) & metaColnames == TRUE & metadata_df == FALSE)) {
    suppressWarnings(metadata_tbl <- readr::read_delim(metadata, col_names = TRUE, delim = delim, col_types = readr::cols()))

  } else if (all(!is.null(metadata) & metaColnames == FALSE & metadata_df == FALSE)) {
    suppressWarnings(metadata_tbl <- readr::read_delim(metadata, col_names = FALSE, delim = delim, col_types = readr::cols()))

    if (length(colnames(metadata_tbl)) != 3) {
      cli::cli_alert_danger('Error: expected three columns: taxonomy, filename, and completeness columns.')
      stop()}

    metadata_tbl <- metadata_tbl %>%
      dplyr::rename_if(~ is.character(.x) & all(!stringr::str_detect(.x, stringr::regex('\\/|\\~|\\.'))), ~ 'taxonomy') %>%
      dplyr::rename_if(~ is.character(.x) & all(stringr::str_detect(.x, stringr::regex('\\/|\\~|\\.'))), ~ 'filename') %>%
      dplyr::rename_if(is.numeric, ~ 'completeness')

  } else if (is.data.frame(metadata_df)) {
    if (all(c('genome_name', 'taxonomy', 'completeness', 'filename') %in% colnames(metadata_df))) {
      metadata_tbl <- metadata_df %>%
        dplyr::select(genome_name, taxonomy, completeness, filename)
    } else {
      cli::cli_alert_danger('Error: Does your metadata dataframe contain the four columns genome_name, taxonomy, completeness, and filename?')
      stop()
    }

  } else {
    stop(cli::cli_alert_danger('Error: metadata file not detected. Please make sure you have provided the correct path to a metadata file.'))
  }

  metadata_tbl <- metadata_tbl %>%
    dplyr::mutate(completeness = completeness / 100)
  if (any(metadata_tbl$completeness > 1 | metadata_tbl$completeness < 0)) {
    cli::cli_alert_danger('Error: Completeness column contains one ore more values with a percentage greater than 100 or less than 0. Please make sure all percentage values are between 0 and 100.')
    stop()
  }

  index <- get_index(.data = Data, metadata_tbl = metadata_tbl)

  if (is.null(Data)) {
    files <- metadata_tbl$filename #list.files(path = filePath, pattern = filePattern, full.names = TRUE)
  }
  if (is.list(Data)) { #colname checks for all forms of Data: NEEDS UPDATE to include all columns that will be REQUIRED
    if (!purrr::every(Data, ~ {
      grepl(paste(colnames(.x), collapse = '|'), c('KO', 'evalue|e-value'), ignore.case = T)})) {
      cli::cli_alert_danger('Error: Colnames KO and e-value were not detected. make sure input columns are properly labelled.')
      stop()
    } else {.}
  }

  formatted_genomes <- list()
  formatted_genomes <- purrr::map(1:index, ~ { #for each element in the index, index is just number vector- # of inputs
    if (is.list(Data) & length(Data) >= 1) { # MAKE SURE THIS ALWAYS RETURNS 1 FOR DATAFRAMES
      Data <- Data[[.x]]
      orgName <- names(Data[.x]) # this is supposed to give a default name each genome tibble in a list if none are provided
    }
    if (!('genome_name' %in% colnames(metadata_tbl))) { #if genome names were not provided...
      if (is.character(Data) | is.null(Data)) {
        files_x <- files[.x]
        #get_fileName.g(filePattern = filePattern, files = files_x)
        orgName <- sub('.*\\/.*\\/+(.*)', '\\1', files_x, perl = TRUE)
        #} else if (is.character(Data)) {
        #  orgName <- get_fileName.g(.data = Data)
      } else if (is.list(Data) & length(Data) >= 1 & !is.data.frame(Data)) {
        if (!is.null(names(Data)[.x])) {
          orgName <- names(Data[.x])
        } else {
          orgName <- paste0('genome_', .x)
        }
      } else if (is.list(Data) & is.data.frame(Data) & length(Data) >= 1) {
        orgName <- deparse(substitute(Data))
      } else {
        cli::cli_alert_danger('Error: Data must be a dataframe, list of dataframes, a flatfile, or multiple flatfiles. Please see usage() for more information.')
        stop()
      }
    } else { # otherwise the above if/else chain is NULL
      orgName <- NULL
    }

    if (is.null(Data)) {
      files_x <- files[.x]
      formatted_genomes[[.x]] <- import_data.g(.data = Data, file = files_x, delim = delim) #%>%
    } else {
      formatted_genomes[[.x]] <- import_data.g(.data = Data, delim = delim)
    }

    formatted_genomes[[.x]] <- formatted_genomes[[.x]] %>%
      {if (kofamscan == TRUE) tidy_kofam(., cutoff = cutoff)
        else if (dram == TRUE) tidy_dram(.)
        else if (kofamscan == FALSE & dram == FALSE) tidy_custom_anno(., cutoff = cutoff, input_type = 'genome_name')
        else stop(cli::cli_alert_danger('Error: Issue importing user data. Please see usage() for data import guidelines.'))}

    if (is.character(Data) | is.null(Data)) {
      metadata_tbl <- metadata_tbl %>%
        dplyr::mutate(filename = stringr::str_replace(filename, '.*\\/.*\\/+(.*)', '\\1')) # this captures the file name w/o full path; does not change after first edit

      formatted_genomes[[.x]] <- formatted_genomes[[.x]] %>%
        dplyr::summarize(k_number = paste0(k_number, collapse = ' '),
                         gene_name = paste0(gene_name, collapse = ' ')) %>%
        {if (!is.null(orgName)) dplyr::mutate(., genome_name = orgName) #genome_name := !!orgName)
          else (.)} %>%
        dplyr::mutate(filename = sub('.*\\/.*\\/+(.*)', '\\1', files_x, perl = TRUE)) %>%
        dplyr::left_join(metadata_tbl, by = 'filename') # make this left_join the metadata row based on a match to the (current) genome file name

    } else {
      formatted_genomes[[.x]] <- formatted_genomes[[.x]] %>%
        concatenate(orgName = orgName, taxon = metadata_tbl$taxonomy[[.x]],
                    completeness = metadata_tbl$completeness[[.x]])
    }
  }) %>%
    dplyr::tibble() %>%
    tidyr::unnest(cols = dplyr::everything()) %>%
    dplyr::group_by(genome_name) %>%
    dplyr::mutate(data_type = 'genome') #%>%
  #dplyr::select(-filename) # should make both methods that have metadata colnames and not colnames use left_join method

  cli::cli_alert_success('Parsed HMM/Blast hits and E-values')
  cli::cli_alert_info('Used E-value cutoff: {cutoff}')

  return(formatted_genomes)
}



get_index <- function(.data, metadata_tbl = metadata_tbl) {
  if (is.character(.data)) {
    index <- length(.data) # could just be put as length 1 right..?
  } else if (is.null(.data)) {
    index <- length(metadata_tbl$filename)
  } else if (is.list(.data) & purrr::every(.data, ~ {grepl(paste(names(.x), collapse = '|'),
                                                           c('KO', `E-value`), ignore.case = TRUE)})) { # NEEDS TO BE CHANGED; KO/K_NUMBER, GENOME_NAME, GENE_NAME, ETC.
    index <- length(.data)
  } else {cli::cli_alert_danger("Error: Issue generating index for data import.")
    stop()
  }
  return(index)
}



import_data.g <- function(.data, file = NULL, delim = NULL) {
  if (is.null(.data)) {
    import <- readr::read_delim(file, col_types = readr::cols(), delim = delim)
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



#get_fileName.g <- function(.data = NULL, filePattern = NULL, files = NULL) {
#  if (is.null(.data)) {fileName <- sub(paste('.*\\/(.*)', filePattern, sep = ''), '\\1', files, perl = T) # files was 'files[.x]'
#  } else if (is.character(.data)) {fileName <- sub('.*\\/(.*)\\..*', '\\1', .data, perl = T) #is.character() may not work as expected here
#  #} else if (is.list(.data)) {fileName <- names(.data)
#  } else {cli::cli_alert_danger('Data input failed. Data must be a dataframe, list of dataframes, or an input flatfile/flatfile(s).')
#    stop()}
#  return(fileName)
#}


#' @export
tidy_kofam <- function(.data, cutoff = 1e-3) {
  .data %>%
    dplyr::filter(!stringr::str_detect(`E-value`, '-----')) %>%
    {if (all(c('E-value', 'KO', 'gene name', 'score') %in% colnames(.))) dplyr::select(., `E-value`, KO, `gene name`, score) #added score
      else stop(cli::cli_alert_danger(
        "Error: Columns 'E-value', 'KO', 'gene name', and/or 'score' not detected. These columns are required to read in Kofamscan output files."))} %>%
    dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
    dplyr::filter(`E-value` <= cutoff | `E-value` == 0) %>%
    dplyr::group_by(`gene name`) %>% ###
    dplyr::filter(score == max(score) & `E-value` == min(`E-value`)) %>% ###
    dplyr::ungroup() %>% ###
    dplyr::select(KO, `gene name`) %>%
    dplyr::rename(gene_name = `gene name`, k_number = KO) %>%
    dplyr::filter(!is.na(k_number), !duplicated(k_number))
}



tidy_dram <- function(.data) {
  .data %>%
    {if (all(c('fasta', 'scaffold', 'gene_position',
               'kegg_id', 'kegg_hit') %in% colnames(.)))  dplyr::select(fasta, scaffold, gene_position, kegg_id, kegg_hit)
      else if (all(c('fasta', 'kegg_id', 'kegg_hit') %in% colnames(.))) dplyr::select(fasta, kegg_id, kegg_hit)
      else if (all(c('kegg_id', 'kegg_hit') %in% colnames(.))) dplyr::select(kegg_id, kegg_hit)
      else stop(cli::cli_alert_danger(
        "Error: Columns 'kegg_id' and 'kegg_hit' not detected. These columns are required to read in Dram output files."))} %>%
    dplyr::filter(!is.na(kegg_id), !duplicated(kegg_id)) %>%
    dplyr::rename(k_number = kegg_id)
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



concatenate <- function(.data, orgName = NULL, taxon = NA, completeness = NA) {
  .data %>%
    #{if (is.null(.data) | is.character(.data)) dplyr::filter(., !(duplicated(k_number))) %>%
    #    dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
    #                     gene_name = paste0(gene_name, collapse = ' ')) %>%
    #    {if (!is.null(orgName)) dplyr::mutate(., genome_name = orgName) #genome_name := !!orgName)
    #      else (.)} %>%
    #    dplyr::mutate(., taxonomy = taxon, completeness = completeness)

    {if (is.data.frame(.data) & length(list(.data)) == 1 & length(.data) >= 1) {
        dplyr::filter(., !(duplicated(k_number))) %>%
          dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
                           gene_name = paste0(gene_name, collapse = ' ')) %>%
          {if (!is.null(orgName)) dplyr::mutate(., genome_name = orgName)  #:= !!deparse(orgName))
            else (.)} %>%
          dplyr::mutate(taxonomy = taxon, completeness = completeness)

      } else {dplyr::filter(., !(duplicated(k_number))) %>% #THIS NEEDS REWORK
          dplyr::summarize(., k_number = paste0(k_number, collapse = ' '),
                           gene_name = paste0(gene_name, collapse = ' ')) %>%
          {if (!is.null(orgName)) dplyr::mutate(., genome_name := !!deparse(paste(orgName, .x, sep = '_'))) # this will fail. needs rework.
            else (.)} %>%
          dplyr::mutate(taxonomy = taxon, completeness = completeness)}
    }
}
