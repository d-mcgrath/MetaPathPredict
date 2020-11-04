#' This function is used to read in user metagenome data gene prediction data with associated taxonomic annotations,
#' and KEGG Orthology HMM/Blast data. Both inputs must be in the form of flatfiles, default is TSV format.
#' @param gene_input The full or relative path to a single flatfile with gene predictions and their associated taxonomic annotations.
#' Columns with taxonomy should include one more of the following: Genus, Family, Order, Class, Phylum, Domain
#' @param ko_input The full or relative path to a  single flatfile with HMM/Blast results which include KEGG Orthology terms in a column titled 'KO'
#' and E-values in a column titled 'E-value'
#' @param evalue The desired E-value cutoff for HMM/Blast KEGG Orthology hits. Default is 0.001
#' @param gene_delim The delimiter of the input gene prediction with taxonomic annotation flatfile. Default is tab
#' @param ko_delim The delimiter of the input KEGG Orthology prediction flatfile. Default is tab
#' @param kofamscan If the input KEGG Orthology flatfile is direct output from Kofamscan, set this argument to TRUE, otherwise FALSE. Default is TRUE
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @export
read_metagenome_data <- function(gene_input, ko_input, evalue = 1e-3, gene_delim = '\t', ko_delim = '\t',
                                 kofamscan = TRUE) {

  taxonomic_lvls <- c('Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain')
  tax_cols <- c('Species' = NA, 'Genus' = NA, 'Family' = NA, 'Order' = NA, 'Class' = NA,
                'Phylum' = NA, 'Domain' = NA)

  cli::cli_h1('Reading metagenomic data into MetaPredict')
  gene_table <- suppressWarnings(readr::read_delim(gene_input, col_types = readr::cols(),
                                                   delim = gene_delim))

  for (level in taxonomic_lvls) {
    gene_table <- gene_table %>%
      rename_at(vars(contains(level)), ~ paste0(level))
  }

  gene_table <- gene_table %>%
    #rename_at(vars(matches(key)), ~ paste0(key)) %>%
    rename_at(vars(contains('Superkingdom')), ~ paste0('Domain')) %>%
    rename_at(vars(contains('ORF')), ~ paste0('ORF')) %>%
    tibble::add_column(!!!tax_cols[!names(tax_cols) %in% names(.)]) %>%
    dplyr::filter(!(Domain == 'not classified'), !(stringr::str_detect(lineage, 'no hit to database')),
                  !(stringr::str_detect(lineage, 'no taxid found')),
                  !(rowSums(across(all_of(taxonomic_lvls), ~ is.na(.))) == 7)) %>%
    mutate(across(all_of(taxonomic_lvls), ~ sub('not classified', NA, .x)),
           organism = coalesce(Genus, Family, Order, Class, Phylum, Domain),
           organism = stringr::str_replace_all(organism,
                                               stringr::regex(c('Candidatus ' = '', '\\*+' = '',
                                                                'candidate division ' = ''),
                                                              ignore_case = TRUE))) %>%
    rename(Gene = ORF) %>%
    select(Gene, organism)

  col <- sub('.*\\/(.*)\\..*', '\\1', ko_input, perl = T)
  ko_table <- suppressWarnings(readr::read_delim(ko_input, col_types = readr::cols(), delim = ko_delim))

  if (kofamscan == TRUE) {
    ko_table <- ko_table %>% slice(-1)
  }
  ko_table <- ko_table %>%
    select(`gene name`, KO, `E-value`) %>%
    mutate(`E-value` = as.numeric(`E-value`)) %>%
    dplyr::filter(`E-value` <= evalue) %>%
    select(`gene name`, KO) %>%
    rename(Gene = `gene name`, ko_term = KO) %>%
    dplyr::filter(!(duplicated(ko_term)))

  gene_table <- gene_table %>%
    full_join(ko_table, by = 'Gene') %>%
    dplyr::filter(!(is.na(ko_term))) %>%
    mutate(organism = case_when(
      is.na(organism) ~ 'unidentified taxonomy',
      TRUE ~ organism)) %>%
    arrange(organism) %>%
    group_by(organism) %>%
    summarize(ko_term = paste0(ko_term, collapse = ' '),
              Gene = paste0(Gene, collapse = ' '),
              organism = unique(organism),
              .groups = 'keep') %>%
    mutate(data_type = 'metagenome',
           metagenome_name := !!col)

  cli::cli_alert_success('Parsed predicted genes and gene taxonomic annotations')
  cli::cli_alert_success('Parsed HMM/Blast hits and E-values')
  cli::cli_alert_info('Used E-value cutoff: {evalue}')

  return(gene_table)
}

