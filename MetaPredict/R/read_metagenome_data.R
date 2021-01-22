#' This function is used to read in user metagenome data gene prediction data with associated taxonomic annotations,
#' and KEGG Orthology HMM/Blast data. Both inputs must be in the form of flatfiles, default is TSV format.
#' @param gene_input The full or relative path to a single flatfile with gene predictions and their associated taxonomic annotations.
#' Columns with taxonomy should include one more of the following: Genus, Family, Order, Class, Phylum, Domain
#' @param ko_input The full or relative path to a  single flatfile with HMM/Blast results which include KEGG Orthology terms in a column titled 'KO'
#' and E-values in a column titled 'E-value'
#' @param covColumn A CSV flatfile containing the estimated metagenome coverage of all input metagenomes
#' @param evalue The desired E-value cutoff for HMM/Blast KEGG Orthology hits. Default is 0.001
#' @param gene_delim The delimiter of the input gene prediction with taxonomic annotation flatfile. Default is tab
#' @param ko_delim The delimiter of the input KEGG Orthology prediction flatfile. Default is tab
#' @param kofamscan If the input KEGG Orthology flatfile is direct output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @importFrom magrittr "%>%"
#' @export
read_metagenome_data <- function(gene_input, ko_input, covColumn = NULL, cutoff = 1e-3, gene_delim = '\t', ko_delim = '\t',
                                 kofamscan = TRUE, cat = TRUE, custom_knumber_ColNames = NULL, custom_anno_ColNames = NULL) {
  if (is.null(covColumn)) {
    stop(cli::cli_alert_danger('Estimated metagenome coverage information required. Please provide a CSV flatfile consisting of a column of the estimated completeness of each metagenome, in the same order as the metagenomes input into MetaPredict. Optionally add a second column with the metagenome name, PROCEEDED by the metagenome coverage (e.g., 1st column: name; 2nd column: coverage).'))
  }
  suppressWarnings(cov_tibble <- readr::read_csv(covColumn,
                                                 col_names = 'coverage', col_types = readr::cols())) #option must be added to have metagenome_name for each metagenome and its coverage

  taxonomic_lvls <- c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')
  tax_cols <- c('Genus' = NA, 'Family' = NA, 'Order' = NA, 'Class' = NA, 'Phylum' = NA, 'Domain' = NA)

  cli::cli_h1('Reading metagenomic data into MetaPredict')
  gene_table <- suppressWarnings(readr::read_delim(gene_input, col_types = readr::cols(),
                                                   delim = gene_delim))
  gene_table <- gene_table %>%
    dplyr::rename_with(~ taxonomic_lvls[stringr::str_detect(taxonomic_lvls, stringr::regex(., ignore_case = TRUE))], .cols = dplyr::contains(taxonomic_lvls)) %>%
    {if (cat == TRUE) tidy_cat_taxa(., taxonomic_lvls = taxonomic_lvls, tax_cols = tax_cols)
      else if (!is.null(custom_anno_ColNames)) tidy_custom_taxa(., cutoff = cutoff)
      else stop(cli::cli_alert_danger('Error: Input format does not match any of the accepted formats. Please see usage()'))}

  ko_table <- suppressWarnings(readr::read_delim(ko_input, col_types = readr::cols(), delim = ko_delim))

  ko_table <- ko_table %>%
    {if (kofamscan == TRUE) tidy_kofam(., cutoff = cutoff, `E-value`, KO, `gene name`)
      else if (!is.null(custom_knumber_ColNames)) tidy_custom_anno(., cutoff = cutoff, input_type = 'metagenome_name')
      else stop(cli::cli_alert_danger('Error message')) }

  if (!('metagenome_name' %in% colnames(cov_tibble)) & !('metagenome_name' %in% colnames(ko_table))) {
    default_name <- sub('.*\\/(.*)\\..*', '\\1', ko_input, perl = T)
    cli::cli_alert_info('Metagenome name(s) not provided/detected. Using default metagenome name: {default_name}')
  } else {default_name <- NULL}

  #need this somewhere .... #dplyr::filter(!(duplicated(k_number)))
  gene_table <- gene_table %>%
    merge_tables(ko_table, cov_tibble, default_name = default_name)

  cli::cli_alert_success('Parsed predicted genes and gene taxonomic annotations')
  cli::cli_alert_success('Parsed HMM/Blast hits and E-values')
  cli::cli_alert_info('Used E-value cutoff: {cutoff}')

  return(gene_table)
}



tidy_cat_taxa <- function(.data, taxonomic_lvls, tax_cols) {
  .data %>%
    dplyr::rename_at(dplyr::vars(dplyr::contains('Superkingdom')), ~ paste0('Domain')) %>%
    dplyr::rename_at(dplyr::vars(dplyr::contains('ORF')), ~ paste0('ORF')) %>%
    tibble::add_column(!!!tax_cols[!names(tax_cols) %in% names(.)]) %>%
    dplyr::mutate(Domain = dplyr::case_when(Domain == 'not classified' ~ NA_character_,
                                            lineage == 'no taxid found' ~ NA_character_,
                                            lineage == 'no hit to database' ~ NA_character_,
                                            TRUE ~ Domain)) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('not classified', NA, .x)),
                  taxonomy = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain), #coalesce takes the lowest tax level as input
                  taxonomy = stringr::str_replace_all(taxonomy,
                                                      stringr::regex(c('Candidatus ' = '', '\\*+' = '',
                                                                       'candidate division ' = ''),
                                                                     ignore_case = TRUE))) %>%
    dplyr::rename(gene_name = ORF) %>%
    dplyr::select(gene_name, taxonomy)
}



#fn to join gene annotation, contig/scaffold taxonomic annotation, and metagenome coverage tibbles
merge_tables <- function(.data, ko_table, cov_tibble, default_name = NULL) {
  .data %>%
    dplyr::full_join(ko_table, by = 'gene_name') %>%
    dplyr::filter(!(is.na(k_number))) %>%
    dplyr::mutate(taxonomy = dplyr::case_when(
      is.na(taxonomy) ~ 'unidentified taxonomy',
      TRUE ~ taxonomy)) %>%
    dplyr::arrange(taxonomy) %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarize(k_number = paste0(k_number, collapse = ' '),
                     gene_name = paste0(gene_name, collapse = ' '),
                     taxonomy = unique(taxonomy),
                     .groups = 'keep') %>%
    {if (!is.null(default_name)) dplyr::mutate(., metagenome_name = default_name)
      else (.)} %>%
    dplyr::mutate(data_type = 'metagenome',# requires if/else statement - names might be provided by user
                  p_j = cov_tibble$coverage[1])
}



tidy_custom_taxa <- function(.data, cutoff = 1e-3) {
  taxonomic_lvls <- c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')
  tax_cols <- c('Genus' = NA, 'Family' = NA, 'Order' = NA, 'Class' = NA, 'Phylum' = NA, 'Domain' = NA)

  .data %>%
    {if ('gene_name' %in% colnames(.) & any(taxonomic_lvls %in% colnames(.))) dplyr::select(., gene_name, dplyr::contains(taxonomic_lvls))
      else stop(cli::cli_alert_danger(
        "Error: Custom taxonomy file must contain column 'gene_name' and at least one column with a taxonomic level (Genus, Family, Order, Class, and/or Phylum).
                The lowest taxonomic level down to Genus will be used. Please see usage() for more information"))} %>%
    purrr::modify(taxonomic_lvls, function(.col) {dplyr::rename(., dplyr::across(dplyr::contains(.col), ~ .col))}) %>%
    tibble::add_column(!!!tax_cols[!names(tax_cols) %in% names(.)]) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('not classified', NA, .x)),
                  taxonomy = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain), #coalesce takes the lowest tax level as input
                  taxonomy = stringr::str_replace_all(taxonomy,
                                                      stringr::regex(c('Candidatus ' = '', '\\*+' = '',
                                                                       'candidate division ' = ''),
                                                                     ignore_case = TRUE))) %>%
    dplyr::select(gene, taxonomy)
}
