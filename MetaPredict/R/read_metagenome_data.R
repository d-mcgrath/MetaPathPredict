#' This function is used to read in user metagenome data gene prediction data with associated taxonomic annotations,
#' and KEGG Orthology HMM/Blast data. Both inputs must be in the form of flatfiles, default is TSV format.
#' @param gene_input The full or relative path to a single flatfile with gene predictions and their associated taxonomic annotations.
#' Columns with taxonomy should include one more of the following: Genus, Family, Order, Class, Phylum, Domain
#' @param ko_input The full or relative path to a  single flatfile with HMM/Blast/other results which include KEGG Orthology terms/K numbers in a column titled 'KO', E-values in a column titled 'E-value', the resulting score is in a column titled 'score', and the gene name is in a column titled 'gene name'
#' @param metadata A flatfile containing the estimated metagenome completeness of the metagenome
#' @param cutoff The desired E-value cutoff for HMM/Blast/other KEGG Orthology/K number hits. Default is 0.001
#' @param gene_df If gene annotations are in a dataframe object, set this argument equal to the name of the object. Default is FALSE
#' @param ko_df If the KEGG Orthology K number annotations are in a dataframe object, set this argument equal to the name of the object. Default is FALSE
#' @param metadata_df If the metadata is in a dataframe object, set this argument equal to the name of the object. Default is FALSE
#' @param gene_delim The delimiter of the input gene prediction with taxonomic annotation flatfile. Default is tab
#' @param ko_delim The delimiter of the input KEGG Orthology prediction flatfile. Default is tab
#' @param metadata_delim The delimiter of the metadata flatfile. Default is tab
#' @param kofamscan If the input KEGG Orthology flatfile is direct output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @param cat If the input gene taxonomic annotations are from CAT, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @importFrom magrittr "%>%"
#' @export
read_metagenome_data <- function(gene_input = NULL, ko_input = NULL, contig_input = NULL, metadata = NULL, gene_df = FALSE, ko_df = FALSE, contig_df = FALSE, metadata_df = FALSE, cutoff = 1e-3, gene_delim = '\t', ko_delim = '\t', contig_delim = '\t', metadata_delim = '\t', kofamscan = TRUE, cat = TRUE, custom_knumber_ColNames = NULL, custom_anno_ColNames = NULL) {

  if (all(!is.null(metadata) & metadata_df == FALSE)) {
    #option must be added to have metagenome_name for each metagenome and its completeness
    suppressWarnings(metadata_tbl <- readr::read_csv(metadata, col_types = readr::cols()))
    if (length(colnames(metadata_tbl)) > 2) {
      cli::cli_alert_danger('Error: expected a maximum of two columns: metagenome_name and completeness columns.')
      stop()
    }
    metadata_tbl <- metadata_tbl %>%
      dplyr::rename_if(~ is.character(.x), ~ 'metagenome_name') %>%
      dplyr::rename_if(is.numeric, ~ 'completeness') %>%
      dplyr::mutate(completeness = completeness / 100)
    if (any(metadata_tbl$completeness > 1 | metadata_tbl$completeness < 0)) {
      cli::cli_alert_danger('Error: completeness column contains one ore more values with a percentage greater than 100 or less than 0. Please make sure all percentage values are between 0 and 100.')
      stop()
    }

  } else if (all(is.null(metadata) & is.data.frame(metadata_df))) {
    if (all(c('metagenome_name', 'completeness') %in% colnames(metadata_df))) {
      metadata_tbl <- metadata_df %>%
        dplyr::select(metagenome_name, completeness) %>%
        dplyr::mutate(completeness = completeness / 100)
      if (any(metadata_tbl$completeness > 1 | metadata_tbl$completeness < 0)) {
        cli::cli_alert_danger('Error: completeness column contains one ore more values with a percentage greater than 100 or less than 0. Please make sure all percentage values are between 0 and 100.')
        stop()
      }
    } else {
      cli::cli_alert_danger('Error: Does your metadata dataframe contain the columns metagenome_name and completeness?')
      stop()
    }

  } else if (all(!is.null(metadata) & is.data.frame(metadata_df))) {
    cli::cli_alert_danger('Error: Either a gene input flatfile or a gene dataframe is required. Both types of inputs were detected. Make sure to use either the metadata or metadata_df argument, not both.')
    stop()

  } else if (all(is.null(metadata) & metadata_df == FALSE)) {
    stop(cli::cli_alert_danger('Metadata including estimated metagenome completeness information required. This should be a value between 0 and 100. Please provide a CSV metadata flatfile or metadata dataframe consisting of a column of the estimated completeness of the metagenome. Optionally add a second column metagenome_name with the metagenome name (column order does not matter).'))

  } else {
    cli::cli_alert_danger('Something went wrong. Please see usage() for instructions to run MetaPredict.')
    stop()
  }

  taxonomic_lvls <- c('Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')
  tax_cols <- c('Genus' = NA, 'Family' = NA, 'Order' = NA, 'Class' = NA, 'Phylum' = NA, 'Domain' = NA)


  # read in gene flatfile, or optionally a dataframe containing the gene information
  if (all(!is.null(gene_input) & gene_df == FALSE)) {
    cli::cli_h1('Reading metagenomic data into MetaPredict')
    gene_table <- suppressWarnings(readr::read_delim(gene_input, col_types = readr::cols(),
                                                     delim = gene_delim))
  } else if (all(!is.null(gene_input) & is.data.frame(gene_df))) {
    cli::cli_alert_danger('Error: Either a gene input flatfile or a gene dataframe is required. Both types of inputs were detected. Make sure to use either the gene_input or gene_df argument, not both.')
    stop()
  } else if (all(is.null(gene_input) & is.data.frame(gene_df))) {
    gene_table <- gene_df
  } else {
    cli::cli_alert_danger('Error: gene input flatfile/gene dataframe not detected. Make sure you have given the path to a flatfile for the gene_input argument, or have provided the name of a dataframe object for the gene_df argument.')
    stop()
  }



  # if using CAT input data, read in contig flatfile, or optionally a dataframe containing the contig information
  if (cat == TRUE) {
    if (all(!is.null(contig_input) & contig_df == FALSE)) {
      cli::cli_h1('Reading contig data into MetaPredict')
      contig_tbl <- suppressWarnings(readr::read_delim(contig_input, col_types = readr::cols(),
                                                       delim = contig_delim))
    } else if (all(!is.null(contig_input) & is.data.frame(contig_df))) {
      cli::cli_alert_danger('Error: Either a contig input flatfile or a contig dataframe is required. Both types of inputs were detected. Make sure to use either the contig_input or contig_df argument, not both.')
      stop()
    } else if (all(is.null(contig_input) & is.data.frame(contig_df))) {
      contig_tbl <- contig_df
    } else {
      cli::cli_alert_danger('Error: contig input flatfile/contig dataframe not detected. Make sure you have given the path to a flatfile for the contig_input argument, or have provided the name of a dataframe object for the contig_df argument.')
      stop()
    }
  }


  gene_table <- gene_table %>%
    dplyr::rename_with(~ taxonomic_lvls[stringr::str_detect(taxonomic_lvls, stringr::regex(., ignore_case = TRUE))], .cols = dplyr::contains(taxonomic_lvls)) %>%
    {if (cat == TRUE) tidy_cat_taxa(.gene_data = ., .contig_data = contig_tbl, taxonomic_lvls = taxonomic_lvls, tax_cols = tax_cols)
      else if (!is.null(custom_anno_ColNames)) tidy_custom_taxa(., cutoff = cutoff)
      else stop(cli::cli_alert_danger('Error: Input format does not match any of the accepted formats. Please see usage()'))}


  # read in gene flatfile, or optionally a dataframe containing the gene information
  if (all(!is.null(ko_input) & ko_df == FALSE)) {
    cli::cli_h1('Reading KEGG Orthology data into MetaPredict')
    ko_table <- suppressWarnings(readr::read_delim(ko_input, col_types = readr::cols(),
                                                   delim = ko_delim))
  } else if (all(!is.null(ko_input) & is.data.frame(ko_df))) {
    cli::cli_alert_danger('Error: Either a KEGG Orthology input flatfile or a KEGG Orthology dataframe is required. Both types of inputs were detected. Make sure to use either one or the other input type, not both.')
    stop()
  } else if (all(is.null(ko_input) & is.data.frame(ko_df))) {
    ko_table <- ko_df
  } else {
    cli::cli_alert_danger('Error: Kegg Orthology input flatfile/gene dataframe not detected. Make sure you have given the path to a flatfile for the ko_input argument, or have provided the name of a dataframe object for the ko_df argument.')
    stop()
  }


  ko_table <- ko_table %>%
    {if (kofamscan == TRUE) tidy_kofam.mg(., cutoff = cutoff)
      else if (!is.null(custom_knumber_ColNames)) tidy_custom_anno(., cutoff = cutoff, input_type = 'metagenome_name')
      else stop(cli::cli_alert_danger('Error message')) }

  if (!('metagenome_name' %in% colnames(metadata_tbl)) & !('metagenome_name' %in% colnames(ko_table))) {
    if (!is.null(ko_input)) {
      default_name <- sub('.*\\/(.*)\\..*', '\\1', ko_input, perl = T)
    } else if (is.data.frame(ko_df)) {
      default_name <- deparse(substitute(ko_df))
    } else {
      stop(cli::cli_alert_danger('Error: Could not create default name for metagenome. Make sure input data and arguments were formatted properly.'))
    }
    cli::cli_alert_info('Metagenome name(s) not provided/detected. Using default metagenome name: {default_name}')
  } else {default_name <- metadata_tbl$metagenome_name}

  #need this somewhere .... #dplyr::filter(!(duplicated(k_number)))
  gene_table <- gene_table %>%
    merge_tables(ko_table, metadata_tbl, default_name = default_name)

  cli::cli_alert_success('Parsed predicted genes and gene taxonomic annotations')
  cli::cli_alert_success('Parsed HMM/Blast hits and E-values')
  cli::cli_alert_info('Used E-value cutoff: {cutoff}')

  return(gene_table)
}



tidy_kofam.mg <- function(.data, cutoff = 1e-3) {
  .data %>%
    dplyr::filter(!stringr::str_detect(`E-value`, '-----')) %>%
    {if (all(c('E-value', 'KO', 'gene name', 'score') %in% colnames(.))) dplyr::select(., `E-value`, KO, `gene name`, score)
      else stop(cli::cli_alert_danger(
        "Error: Columns 'E-value', 'KO', 'gene name', and/or 'score' not detected. These columns are required to read in Kofamscan output files."))} %>%
    dplyr::mutate(`E-value` = as.numeric(`E-value`)) %>%
    dplyr::filter(`E-value` <= cutoff | `E-value` == 0) %>%
    dplyr::group_by(`gene name`) %>%
    dplyr::filter(score == max(score) & `E-value` == min(`E-value`)) %>%
    dplyr::ungroup() %>%
    dplyr::select(KO, `gene name`) %>%
    dplyr::rename(gene_name = `gene name`, k_number = KO) %>%
    dplyr::filter(!is.na(k_number))
}



#tidy_cat_taxa <- function(.data, taxonomic_lvls, tax_cols) {
#  .data %>%
#    dplyr::rename_at(dplyr::vars(dplyr::contains('Superkingdom')), ~ paste0('Domain')) %>%
#    dplyr::rename_at(dplyr::vars(dplyr::contains('ORF')), ~ paste0('gene_name')) %>%
#    tibble::add_column(!!!tax_cols[!names(tax_cols) %in% names(.)]) %>%
#    dplyr::mutate(Domain = dplyr::case_when(Domain == 'not classified' ~ NA_character_,
#                                            lineage == 'no taxid found' ~ NA_character_,
#                                            lineage == 'no hit to database' ~ NA_character_,
#                                            TRUE ~ Domain)) %>%
#    dplyr::mutate(dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('not classified', NA, .x)),
#                  dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('.*\\*$', NA, .x)),
#                  taxonomy = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain)) %>% #, #coalesce takes the lowest tax level as input
    #taxonomy = stringr::str_replace_all(taxonomy,
    #                                    stringr::regex(c('Candidatus ' = '', '\\*+' = '',
    #                                                     'candidate division ' = ''),
    #                                                   ignore_case = TRUE))) %>%
    #dplyr::rename(gene_name = ORF) %>%
#    dplyr::select(gene_name, taxonomy) %>%
#    dplyr::filter(taxonomy != 'Viruses')
#}



# works, but may need to add function to filter out any genes that dont match to contigs
tidy_cat_taxa <- function(.gene_data, .contig_data, taxonomic_lvls, tax_cols) {
  .gene_data %>%
    dplyr::rename_at(dplyr::vars(dplyr::contains('ORF')), ~ paste0('gene_name')) %>%
    dplyr::select(gene_name) %>%
    dplyr::mutate(contig = str_replace(gene_name, '(.*)_\\d+', '\\1')) %>%
    dplyr::left_join(.contig_data %>%
                       dplyr::rename_at(dplyr::vars(dplyr::contains('contig')), ~ paste0('contig')) %>%
                       dplyr::rename_at(dplyr::vars(dplyr::contains('Superkingdom')), ~ paste0('Domain')) %>%
                       dplyr::rename_at(dplyr::vars(dplyr::contains(taxonomic_lvls)), ~ str_to_title(.x)) %>%
                       tibble::add_column(!!!tax_cols[!names(tax_cols) %in% names(.)]) %>%
                       dplyr::mutate(dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('not classified', NA, .x)),
                                     dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('no taxid found', NA, .x)),
                                     dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('no hit to database', NA, .x)),
                                     dplyr::across(dplyr::all_of(taxonomic_lvls), ~ stringr::str_replace(.x, '(.*):.*', '\\1')),
                                     dplyr::across(dplyr::all_of(taxonomic_lvls), ~ sub('.*\\*$', NA, .x))),
                     by = 'contig') %>%
    dplyr::mutate(taxonomy = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain)) %>% #coalesce takes the lowest tax level as input
    dplyr::select(gene_name, taxonomy) %>%
    dplyr::filter(!(taxonomy %in% 'Eukaryota'),
                  !(taxonomy %in% 'Viruses')) #filter out genes from contigs classified as viral or eukaryotic
}



#fn to join gene annotation, contig/scaffold taxonomic annotation, and metagenome completeness tibbles
merge_tables <- function(.data, ko_table, metadata_tbl, default_name = NULL) {
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
      else stop(cli::cli_alert_danger('Error: Metagenome name not detected. Did you include a metagenome_name column in your metadata flatfile/dataframe?'))} %>%
    dplyr::mutate(data_type = 'metagenome',# requires if/else statement - names might be provided by user
                  completeness = metadata_tbl$completeness[1])
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
                  taxonomy = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain)) %>% #, #coalesce takes the lowest tax level as input
    #taxonomy = stringr::str_replace_all(taxonomy,
    #                                    stringr::regex(c('Candidatus ' = '', '\\*+' = '',
    #                                                     'candidate division ' = ''),
    #                                                   ignore_case = TRUE))) %>%
    dplyr::select(gene, taxonomy)
}
