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

  #add something like if (gene_input %>% contains(any_of(taxonomic levels))) {}...
  message('Parsing predicted genes and gene taxonomic annotations into MetaPredict')
  gene_table <- suppressWarnings(readr::read_delim(gene_input, col_types = readr::cols(),
                                  delim = gene_delim)) %>%
    dplyr::filter(!(is.na(superkingdom) | superkingdom == 'not classified'),
           !(stringr::str_detect(lineage, 'no taxid found'))) %>%
    rename(Gene = `# ORF`, Phylum = phylum, Class = class, Order = order, Family = family,
           Genus = genus, Species = species, Domain = superkingdom) %>% #this renaming will need to be made conditional to accept input from other tools besides CAT
    mutate(across(.cols = c(Genus, Family, Order, Class, Phylum, Species), ~ sub('not classified', NA, .x)),
           organism = coalesce(Genus, Family, Order, Class, Phylum, Domain),
           organism =
             stringr::str_replace_all(organism,
                                      stringr::regex(c('Candidatus ' = '', '\\*+' = '',
                                                       'candidate division ' = ''),
                                                     ignore_case = TRUE))) %>%
    select(Gene, organism)

  message('Parsing HMM/Blast hits and E-values into MetaPredict. Using E-value cutoff: ', evalue)
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
    arrange(organism)

  present_na <- gene_table %>% dplyr::filter(is.na(organism))

  gene_table <- gene_table %>%
    dplyr::filter(!is.na(organism)) %>%
    group_by(organism)

  res <- list('gene_table' = gene_table, 'present_na' = present_na)
  return(res)
}

