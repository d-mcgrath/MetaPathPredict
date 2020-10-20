#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param userData Created with the read_data function - KEGG Orthology data for one or more bacteria/archaea
#' @param taxonList A CSV flatfile with lowest taxonomic level of each input organism. Same order as in userData, one per line
#' @param cores Number of cores to use for computation. Default is 1
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @export
MetaPredict <- function(userData, taxonList = NULL, cores = 1) {
  if (!is.null(taxonList)) {
    suppressWarnings(taxa <- readr::read_csv(taxonList,
                                             col_names = 'col', col_types = readr::cols()))
  }
  if (is.null(taxonList)) {
    na_add <- userData[[2]] %>%
      inner_join(ko_term.tibble, by = 'ko_term') %>%
      dplyr::filter(!(duplicated(reaction))) %>%
      inner_join(pathways.tibble, by = c('reaction', 'reaction_description')) %>%
      arrange(pathway) %>%
      mutate(organism = 'unidentified taxonomy', probability = 'Present')

    userData <- userData[[1]]
  }
  temp <- userData %>%
    inner_join(ko_term.tibble, by = 'ko_term') %>%
    dplyr::filter(!(duplicated(reaction))) #at this point we have ko_term, reaction, and reaction_description columns

  res <- temp %>%
    group_split() %>%
    purrr::map(full_join, pathways.tibble, temp, by = c('reaction', 'reaction_description')) %>%
    purrr::map(arrange, pathway) %>%
    purrr::map(group_by, pathway) %>%
    purrr::map(dplyr::filter, !(all(is.na(organism))), !(is.na(pathway))) %>%
    purrr::keep(~ !all(is.na(.x$reaction))) #filters out empty lists

  scan_missing <- res %>%
    purrr::map(ungroup) %>%
    purrr::map(dplyr::filter, is.na(ko_term), reaction %in% colnames(bacteria.rxn.matrix)) %>% #bacteria/archaea matrices have the same column names
    purrr::map(select, -ko_term)

  message('Setting parallel computing parameters, using ', cores, ' core(s)...')
  mplan <- future::plan('multicore', workers = cores) # it is recommended not to call plan() in a function... leave it up to the user?..
  on.exit(future::plan(mplan), add = TRUE)
  options(future.globals.maxSize = 3145728000)

  if (!is.null(taxonList)) {
    results <- purrr::map(1:length(taxa$col), ~ { #added scan here...changed 'predictions' to 'results'..
      message('Processing ', unique(userData$organism)[.x], '...')
      pull_data(scan_missing[[.x]]$reaction, taxa$col[[.x]], scan_missing[[.x]], # NEED to ADJUST THIS to allow for input from metagenome function
                res[[.x]]) }) # added .y (scan) here...

  } else {
    results <- purrr::map(1:length(scan_missing), ~ { #may want to use furrr::future_map() here for optional parallel processing
      message('Processing ', unique(na.omit(res[[.x]]$organism)), ' gene annotations...')
      if (length(scan_missing[[.x]]$reaction) > 0) {
        pull_data(scan_missing[[.x]]$reaction, unique(na.omit(res[[.x]]$organism)),
                  scan_missing[[.x]], # NEED to ADJUST THIS to allow for input from metagenome function
                  res[[.x]])
      } else {
        message('Skipping ', unique(res[[.x]]$organism), '. Nothing to scan...')
        return(res[[.x]]) #return just res[[.x]] since there is nothing to scan...
        }
      })
  }
  if (is.null(taxonList)) {
    results <- results %>%
      bind_rows() %>%
      mutate(probability = case_when(
        !is.na(reaction) & is.na(probability) ~ 'reaction has no associated KEGG Orthology term(s); no calculation done',
        TRUE ~ probability)) %>%
      bind_rows(na_add) %>%
      arrange(pathway, pathway_step)
    }
  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')
  message('Writing results and summary file. Output is in a list. See [[1]] for the summary, [[2]] for the full results.')

  summary_information <- results %>%
    dplyr::filter(!(duplicated(reaction))) %>%
    tally(probability == 'Present') %>%
    full_join(results %>%
                group_by(pathway, reaction) %>%
                dplyr::filter(!any(probability == 'Present')) %>%
                mutate(probability = as.numeric(probability)) %>%
                dplyr::filter(!is.na(probability)) %>%
                ungroup() %>%
                group_by(pathway) %>%
                dplyr::filter(!(duplicated(reaction))) %>%
                tally(length(pathway)), by = 'pathway') %>%
    full_join(results %>%
                dplyr::filter(!(duplicated(reaction))) %>%
                tally(length(pathway)), by = 'pathway') %>%
    rename(pathway_steps_present = 2, pathway_steps_predicted = 3,
           total_pathways_steps = 4) %>%
    mutate(pathway_steps_predicted = case_when(
      is.na(pathway_steps_predicted) ~ 0L,
      TRUE ~ pathway_steps_predicted)) %>%
    inner_join(select(results, pathway, pathway_name), by = 'pathway') %>%
    select(pathway, pathway_name, pathway_steps_present,
           pathway_steps_predicted, total_pathways_steps)

  return(list(summary_information, results))
  }
