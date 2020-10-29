#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param userData Created with the read_data function - KEGG Orthology data for one or more bacteria/archaea
#' @param taxonList A CSV flatfile with lowest taxonomic level of each input organism. Same order as in userData, one per line
#' @param cores Number of cores to use for computation. Default is 1
#' @param output_dir The full or relative path to an output directory where result and summary output will be saved as TSV flatfiles
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @export
new_MetaPredict <- function(userData, taxonList = NULL, cores = 1, output_dir = NULL) {
  if (!is.null(taxonList)) {
    suppressWarnings(taxa <- readr::read_csv(taxonList,
                                             col_names = 'col', col_types = readr::cols()))
  }

  message('Formatting input data...')
  tidy_input <- format_input(userData)
  present_modules <- pull_present_modules_from(tidy_input, userData)

  scan_missing <- present_modules %>%
    purrr::map(ungroup) %>%
    purrr::map(dplyr::filter, is.na(probability))

  #add conditional here to only execute these commands if parallel computing setting is TRUE
  #message('Setting parallel computing parameters, using ', cores, ' core(s)...')
  #mplan <- future::plan('multicore', workers = cores) # it is recommended not to call plan() in a function... leave it up to the user?..
  #on.exit(future::plan(mplan), add = TRUE)
  #options(future.globals.maxSize = 3145728000)


  #change the code below this line to make the updated function work properly
  if (!is.null(taxonList)) {
    results <- purrr::map(1:length(taxa$col), ~ {
      message('Processing ', unique(userData$organism)[.x], '...')
      pull_data(scan_missing[[.x]]$step, taxa$col[[.x]], scan_missing[[.x]],
                present_modules[[.x]])
    })

  } else {
    results <- purrr::map(1:length(scan_missing), ~ { #may want to use furrr::future_map() here for optional parallel processing
      message('Processing ', unique(present_modules[[.x]]$organism), ' gene annotations...')

      if (length(scan_missing[[.x]]$step) > 0 & userData$organism[.x] != 'unidentified taxonomy') {
        pull_data(scan_missing[[.x]]$step, unique(present_modules[[.x]]$organism),
                  scan_missing[[.x]],
                  present_modules[[.x]])

      } else if (length(scan_missing[[.x]]$step) == 0 & userData$organism[.x] != 'unidentified taxonomy') {
        message('Skipping ', userData$organism[.x], '. Nothing to scan...')
        return(present_modules[[.x]]) #return just present_modules[[.x]] since there is nothing to scan...

      } else {
        message('Skipping ', userData$organism[.x], '. Calculations for unidentified taxonomies currently not supported...')
        }
    })
  }
  if (is.null(taxonList)) {
    results <- results %>%
      bind_rows() %>%
      arrange(module, step)
  }
  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')
  message('Writing results and summary file. Output is in a list. Enter "View(results[[{x}]][[1]]) to look at summary information, View(results[[{x}]][[2]] for the full results, and View(results[[{x}]][[3]][[{y}]]) for a heatmaps; change {x} to a number of one or more single genomes to view results for a particular genome. Change {y} to view heatmaps; e.g., enter 1 to view the first heatmap. For single genome(s), pathways are sorted by their size from the smallest in the first heatmap to the largest in the final heatmap.')

  if (is.null(taxonList)) {
    summary_information <- summarize_output(results)
    return(list(summary_information, results))

  } else {
    summary_information <- purrr::map(results, summarize_output)
    heatmaps <- purrr::map(results, create_heatmaps_from)

    return(purrr::transpose(list(summary_information, results, heatmaps)))
  }
}



format_input <- function(userData) {
  input_data <- userData %>%
    select(-Gene) %>%
    tibble::column_to_rownames(var = 'organism') %>%
    t() %>%
    as_tibble()

  input_data <- input_data[rep(seq_len(nrow(input_data)), each = length(nrow(patt.kegg_modules))), ]

  input_formatted <- patt.kegg_modules %>%
    bind_cols(input_data) %>%
    group_by(step, patt_num) %>%
    summarize(across(everything(), ~ all(stringr::str_detect(
      string = .x, pattern = patt))), .groups = 'keep') %>%
    ungroup(patt_num) %>%
    select(-c(name, module, patt, patt_num)) %>%
    summarize(across(everything(), ~ any(.x)), .groups = 'drop') %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'organism') %>%
    slice(-1) %>%
    group_by(organism) %>%
    group_split() %>%
    purrr::map(~ {
      .x <- .x %>%
        t() %>%
        as_tibble() %>%
        slice(-1)
    })
  return(input_formatted)
}



pull_present_modules_from <- function(tidy_data, userData) {
  present_modules <- list()
  present_modules <- purrr::map(1:length(tidy_data), ~ {
    present_modules[[.x]] <- concat.kegg_modules %>%
      ungroup() %>%
      mutate(to_keep = tidy_data[[.x]],
             probability = 'Present') %>%
      dplyr::filter(to_keep == TRUE) %>%
      select(-to_keep) %>%
      full_join(concat.kegg_modules, by = c('name', 'module', 'full', 'step')) %>%
      tibble::add_column(organism = userData$organism[.x], .before = 1) %>%
      arrange(step) %>%
      group_by(module) %>%
      dplyr::filter(!all(is.na(probability))) # filters out pathways for which there is no evidence of their presence in the gene annotations;
    # this should ultimately be a setting that can be turned on/off...
  })
  return(present_modules)
}



get_splits_from <- function(heatmap_matrices, n_chunks) {
  keys <- split(unique(heatmap_matrices$module),
                sort(1:length(unique(heatmap_matrices$module)) %% n_chunks))
  return(keys)
}



determine_splits_for <- function(heatmap_matrices) {
  if (n_groups(heatmap_matrices) <= 77) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 1)

  } else if (n_groups(heatmap_matrices) > 77 & n_groups(heatmap_matrices) <= 154) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 2)

  } else if (n_groups(heatmap_matrices) > 154 & n_groups(heatmap_matrices) <= 231) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 3)

  } else if (n_groups(heatmap_matrices) > 231 & n_groups(heatmap_matrices) <= 308) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 4)

  } else if (n_groups(heatmap_matrices) > 308 & n_groups(heatmap_matrices) <= 385) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 5)

  } else if (n_groups(heatmap_matrices) > 385 & n_groups(heatmap_matrices) <= 462) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 6)

  } else if (n_groups(heatmap_matrices) > 462 & n_groups(heatmap_matrices) <= 539) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 7)

  } else if (n_groups(heatmap_matrices) > 539 & n_groups(heatmap_matrices) <= 616) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 8)

  } else if (n_groups(heatmap_matrices) > 616 & n_groups(heatmap_matrices) <= 693) {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 9)

  } else {
    splits <- get_splits_from(heatmap_matrices, n_chunks = 10)
  }
  return(splits)
}



create_heatmaps_from <- function(results) {
  heatmap_matrices <- results %>%
    mutate(probability = case_when(probability == 'Present' ~ '1',
                                   TRUE ~ probability),
           step = paste0(seq_along(module)),
           module_length = length(module)) %>%
    arrange(module_length) %>%
    select(-c(full, module_length))

  splits <- determine_splits_for(heatmap_matrices)

  heatmap_matrices <- heatmap_matrices %>%
    ungroup() %>%
    mutate(split = purrr::map(1:length(splits), ~ case_when(module %in% splits[[.x]] ~ .x)) %>%
             purrr::as_vector() %>%
             na.omit()) %>%
    group_by(split) %>%
    group_split() %>%
    purrr::map(~ {
      .x <- .x %>%
        ungroup() %>%
        select(name, step, probability) %>%
        tidyr::pivot_wider(id_cols = step, names_from = name, values_from = probability,
                    values_fn = toString) %>%
        select(-step) %>%
        mutate(across(everything(), as.numeric)) %>%
        as.matrix() %>%
        t()
    })

  heatmap_list <- purrr::map(heatmap_matrices, ~ {
    pheatmap::pheatmap(.x, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE,
                       cellwidth = 10, cellheight = 9, color = colorRampPalette(c('snow', 'seagreen4'))(100))
  })
  return(heatmap_list)
}



#  if (!is.null(output_dir) & is.null(taxonList)) {
#    purrr::map(list(summary_information, results), ~ {
#      readr::write_tsv(.x, file = )
#    })
#  }

#  if (!is.null(output_dir) & !is.null(taxonList)) {
#    purrr::map(purrr::transpose(list(summary_information, results)), ~ {
#      readr::write_tsv(.x[[1]][[1]], file = )
#    })
#  }
#}

#purrr::map2(1:length(results[[1]]), c('summary', 'results'), ~ {
#  readr::write_tsv(results[[1]][[.x]], paste0('~/Downloads/metapredict-', .y, '.tsv'))
#})


