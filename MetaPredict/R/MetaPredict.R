#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#'
#' @param userData Created with the read_data function - KEGG Orthology data for one or more bacteria/archaea
#' @param taxonList A CSV flatfile with lowest taxonomic level of each input organism. Same order as in userData, one per line
#' @param cores Number of cores to use for computation. Default is 1
#' @importFrom magrittr "%>%"
#' @export
MetaPredict <- function(userData, taxonList, cores = 1) {
  orgs <- readr::read_csv(taxonList, col_names = 'col', col_types = readr::cols())

  temp <- userData %>%
    dplyr::inner_join(ko_term.tibble, by = 'ko_term') %>%
    dplyr::filter(!(duplicated(reaction))) #at this point we have ko_term, reaction, and reaction_description columns

  res <- temp %>%
    dplyr::group_split() %>%
    purrr::map(dplyr::full_join, pathways.tibble, temp, by = c('reaction', 'reaction_description')) %>%
    purrr::map(dplyr::arrange, pathway) %>%
    purrr::map(dplyr::group_by, pathway) %>%
    purrr::map(dplyr::filter, !(all(is.na(organism))), !(is.na(pathway)))

  scan_missing <- res %>%
    purrr::map(dplyr::ungroup) %>%
    purrr::map(dplyr::filter, is.na(ko_term), reaction %in% colnames(bacteria.rxn.matrix)) %>% #bacteria/archaea matrices have the same column names
    purrr::map(dplyr::select, -ko_term)

  message('Setting parallel computing parameters, using ', cores, ' core(s)...')
  mplan <- future::plan('multicore', workers = cores) # it is recommended not to call plan() in a function... leave it up to the user?..
  on.exit(future::plan(mplan), add = TRUE)
  options(future.globals.maxSize = 3145728000)

  results <- purrr::map(1:length(orgs$col), ~ { #added scan here...changed 'predictions' to 'results'..
    message('Processing ', unique(userData$organism)[.x], '...')
    pull_data(scan_missing[[.x]]$reaction, orgs$col[[.x]], scan_missing[[.x]],
              res[[.x]]) }) # added .y (scan) here...

  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')

  return(results)
}
