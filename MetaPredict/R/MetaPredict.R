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
    results <- purrr::map(1:length(taxa$col), ~ {
      message('Processing ', unique(userData$organism)[.x], '...')
      pull_data(scan_missing[[.x]]$reaction, taxa$col[[.x]], scan_missing[[.x]],
                res[[.x]])
        }) %>%
      purrr::map(~ {
        .x <- .x %>%
          mutate(probability = case_when(
            !is.na(reaction) & is.na(probability) ~ 'Not present; reaction has no associated KEGG Orthology term(s); no calculation done',
            TRUE ~ probability)) %>%
          arrange(pathway, pathway_step)
      })


  } else {
    results <- purrr::map(1:length(scan_missing), ~ { #may want to use furrr::future_map() here for optional parallel processing
      message('Processing ', unique(na.omit(res[[.x]]$organism)), ' gene annotations...')

      if (length(scan_missing[[.x]]$reaction) > 0) {
        pull_data(scan_missing[[.x]]$reaction, unique(na.omit(res[[.x]]$organism)),
                  scan_missing[[.x]],
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
        !is.na(reaction) & is.na(probability) ~ 'Not present; reaction has no associated KEGG Orthology term(s); no calculation done',
        TRUE ~ probability)) %>%
      bind_rows(na_add) %>%
      arrange(pathway, pathway_step)
    }
  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')
  message('Writing results and summary file. Output is in a list. See [[1]] for the summary, [[2]] for the full results.')

  if (is.null(taxonList)) {
    summary_information <- summarize_output(results)
    return(list(summary_information, results))

  } else {
      summary_information <- purrr::map(results, summarize_output)
      return(purrr::transpose(list(summary_information, results)))
    }
  }
