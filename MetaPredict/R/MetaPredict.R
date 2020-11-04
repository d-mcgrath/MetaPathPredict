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
MetaPredict <- function(userData, taxonList = NULL, cores = 1, output_dir = NULL) {
  if (!is.null(taxonList)) {
    suppressWarnings(taxa <- readr::read_csv(taxonList,
                                             col_names = 'col', col_types = readr::cols()))
  }

  cli::cli_h1('Starting MetaPredict')
  cli::cli_alert_info('Formatting input data...')

  if (missing(userData)) {
    cli::cli_alert_danger('Input object not detected in global environment. Make sure you have run read_genome_data() or read_metagenome_data() on your data.')
    stop()
  }
  if (unique(userData$data_type) == 'metagenome' & !is.null(taxonList)) {
    cli::cli_alert_danger('Input is from a metagenome. A taxon list is not required for this data type. Try running MetaPredict() without the taxonList argument.')
    stop()
  }
  if (unique(userData$data_type) == 'genome' & is.null(taxonList)) {
    cli::cli_alert_danger('Input is from one or more genomes. A taxon list is required for this data type. Try running MetaPredict() with the taxonList argument.')
    stop()
  }

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

  cli::cli_alert_info('Reconstructing metabolic pathways and performing prediction calculations...')

  #change the code below this line to make the updated function work properly
  if (!is.null(taxonList)) {
    results <- purrr::map(1:length(taxa$col), ~ {
      #cli::cli_alert_info('Processing {unique(userData$organism)[.x]} ...')
      pull_data(scan_missing[[.x]]$step, taxa$col[[.x]], scan_missing[[.x]],
                present_modules[[.x]])
    }) %>%
      purrr::set_names(unique(userData$organism))

  } else {
    results <- purrr::map(1:length(scan_missing), ~ { #may want to use furrr::future_map() here for optional parallel processing
      if (length(scan_missing[[.x]]$step) > 0 & userData$organism[.x] != 'unidentified taxonomy') {
        pull_data(scan_missing[[.x]]$step, unique(present_modules[[.x]]$organism),
                  scan_missing[[.x]],
                  present_modules[[.x]])
        #cli::cli_alert_success('Processed {unique(present_modules[[.x]]$organism)} gene annotations')

      } else if (length(scan_missing[[.x]]$step) == 0 & userData$organism[.x] != 'unidentified taxonomy') {
        #cli::cli_alert_info('Skipping {userData$organism[.x]}. Nothing to scan')
        return(present_modules[[.x]]) #return just present_modules[[.x]] since there is nothing to scan...

      } else {
        #cli::cli_alert_warning('Skipping {userData$organism[.x]}. Calculations for unidentified taxonomies currently not supported')
        }
    })
  }
  if (is.null(taxonList)) {
    results <- results %>%
      bind_rows() %>%
      arrange(module, step) %>%
      mutate(metagenome_name = unique(userData$metagenome_name))
  }
  if (is.null(taxonList)) {
    summary_information <- summarize_output(results)
    heatmaps <- create_heatmaps_from(results, userData, metagenome = TRUE)
    output <- list(list(summary_information, results, heatmaps)) %>%
      purrr::set_names(unique(userData$metagenome_name))

    if (is.null(output_dir)) {
      cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and reaction probability calculations. Output is in a list.')
      cli::cli_alert_info('Enter View(results[[(x)]][[1]]) to look at summary information, View(results[[(x)]][[2]] for the full results, and print(results[[(x)]][[3]]) for pathway completion heatmaps.')
      cli::cli_alert_warning('Make sure to change (x) to a number between 1 and x, with x equal to the total amount of input genome annotations.')
      cli::cli_alert_info("To save results, run the following command: save_results(output_dir = '/path/to/output/directory'). Note: the output directory will be created if it does not exist.")
      return(output)
    }
  } else {
    summary_information <- purrr::map(results, summarize_output)
    heatmaps <- purrr::map(results, ~ {create_heatmaps_from(.x, userData)})
    output <- purrr::transpose(list(summary_information, results, heatmaps))

    if (is.null(output_dir)) {
      cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and reaction probability calculations. Output is in a list.')
      cli::cli_alert_info('Enter View(results[[(x)]][[1]]) to look at summary information, View(results[[(x)]][[2]] for the full results, and print(results[[(x)]][[3]]) for pathway completion heatmaps.')
      cli::cli_alert_warning('Make sure to change (x) to a number between 1 and x, with x equal to the total amount of input genome annotations.')
      cli::cli_alert_info("To save results, run the following command: save_results(output_dir = '/path/to/output/directory/'). Note: the output directory will be created if it does not exist.")
      return(output)
    }
  }
  if (!is.null(output_dir)) {
    save_results(output, output_dir)
    cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and reaction probability calculations. Output is in directory: {output_dir}')
    return(output)
  }
}



format_input <- function(userData) {
  input_data <- userData
  if ('metagenome_name' %in% colnames(userData)) {
    input_data <- userData %>% select(-metagenome_name)
  }
  input_data <- input_data %>%
    select(-c(Gene, data_type)) %>%
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
      dplyr::filter(!all(is.na(probability))) # filters out pathways for which there is no evidence of their presence in the gene annotations; this should ultimately be a setting that can be turned on/off...
    #if ('metagenome_name' %in% colnames(userData)) {
    #  present_modules[[.x]] <- present_modules[[.x]] %>%
    #    mutate(metagenome_name = unique(userData$metagenome_name))
    #}
  })
  return(present_modules)
}



get_splits_from <- function(heatmaps, n_chunks) {
  keys <- split(unique(heatmaps$module),
                ggplot2::cut_number(1:length(unique(heatmaps$module)), n_chunks)) %>%
    purrr::set_names(1:n_chunks)
  return(keys)
}



determine_splits_for <- function(heatmaps, int = 51, n_ints = 1:15) {
  ints = int * n_ints

  if (n_groups(heatmaps) <= ints[1]) {
    splits <- get_splits_from(heatmaps, n_chunks = 1)

  } else if (n_groups(heatmaps) > ints[1] & n_groups(heatmaps) <= ints[2]) {
    splits <- get_splits_from(heatmaps, n_chunks = 2)

  } else if (n_groups(heatmaps) > ints[2] & n_groups(heatmaps) <= ints[3]) {
    splits <- get_splits_from(heatmaps, n_chunks = 3)

  } else if (n_groups(heatmaps) > ints[3] & n_groups(heatmaps) <= ints[4]) {
    splits <- get_splits_from(heatmaps, n_chunks = 4)

  } else if (n_groups(heatmaps) > ints[4] & n_groups(heatmaps) <= ints[5]) {
    splits <- get_splits_from(heatmaps, n_chunks = 5)

  } else if (n_groups(heatmaps) > ints[5] & n_groups(heatmaps) <= ints[6]) {
    splits <- get_splits_from(heatmaps, n_chunks = 6)

  } else if (n_groups(heatmaps) > ints[6] & n_groups(heatmaps) <= ints[7]) {
    splits <- get_splits_from(heatmaps, n_chunks = 7)

  } else if (n_groups(heatmaps) > ints[7] & n_groups(heatmaps) <= ints[8]) {
    splits <- get_splits_from(heatmaps, n_chunks = 8)

  } else if (n_groups(heatmaps) > ints[8] & n_groups(heatmaps) <= ints[9]) {
    splits <- get_splits_from(heatmaps, n_chunks = 9)

  } else if (n_groups(heatmaps) > ints[9] & n_groups(heatmaps) <= ints[10]) {
    splits <- get_splits_from(heatmaps, n_chunks = 10)

  } else if (n_groups(heatmaps) > ints[10] & n_groups(heatmaps) <= ints[11]) {
    splits <- get_splits_from(heatmaps, n_chunks = 11)

  } else if (n_groups(heatmaps) > ints[11] & n_groups(heatmaps) <= ints[12]) {
    splits <- get_splits_from(heatmaps, n_chunks = 12)

  } else if (n_groups(heatmaps) > ints[12] & n_groups(heatmaps) <= ints[13]) {
    splits <- get_splits_from(heatmaps, n_chunks = 13)

  } else if (n_groups(heatmaps) > ints[13] & n_groups(heatmaps) <= ints[14]) {
    splits <- get_splits_from(heatmaps, n_chunks = 14)

  } else {
    splits <- get_splits_from(heatmaps, n_chunks = 15)
  }
  return(splits)
}



create_heatmaps_from <- function(results, userData, metagenome = FALSE) {
  if (metagenome == FALSE) {
    heatmaps <- results %>%
      mutate(probability = case_when(probability == 'Present' ~ '1',
                                     TRUE ~ probability),
             probability = as.numeric(probability),
             step = stringr::str_replace(step, 'M\\d{5}_', ''),
             module_length = length(module)) %>%
      arrange(module_length) %>%
      select(-c(full, module_length))

    title_org <- unique(userData$organism)

  } else {
    heatmaps <- results %>%
      filter(probability != 'taxonomy not found; no calculation done') %>%
      mutate(probability = case_when(probability == 'Present' ~ '1',
                                     TRUE ~ probability),
             probability = as.numeric(probability)) %>%
      group_by(module, name, step) %>%
      summarize(probability = max(probability), .groups = 'drop') %>%
      group_by(module) %>%
      mutate(step = stringr::str_replace(step, 'M\\d{5}_', ''),
             module_length = length(module)) %>%
      arrange(module_length) %>%
      select(-module_length)

    title_org <- unique(userData$metagenome_name)
  }
  splits <- determine_splits_for(heatmaps)

  heatmaps <- heatmaps %>%
    ungroup() %>%
    mutate(split = purrr::map(1:length(splits), ~ case_when(module %in% splits[[.x]] ~ .x)) %>%
             purrr::as_vector() %>%
             na.omit()) %>%
    mutate(step = factor(step, levels = unique(step)), name = factor(name, levels = unique(name))) %>%
    ggplot2::ggplot(ggplot2::aes(step, forcats::fct_rev(name), fill = probability)) +
    ggplot2::labs(x = 'Module step', y = 'Module name', fill = 'Pathway completion probability',
                  title = paste0(title_org, ': Metabolic pathway predictions')) +
    ggplot2::geom_tile(color = 'snow') +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(low = 'lightskyblue', mid = '#FFFF99', high = 'lightcoral', midpoint = 0.5) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4), axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~ split, scales = 'free')

  return(heatmaps)
}



#' @export
save_results <- function(output, output_dir) {
  if (!dir.exists(output_dir)) {
    cli::cli_alert_warning('Creating output directory {output_dir}')
    dir.create(output_dir)
  }
  if (!stringr::str_detect(output_dir, '.*/$')) {
    output_dir <- sub('(.*)', '\\1\\/', output_dir, perl = TRUE)
  }
  purrr::map(1:length(output), ~ {
    readr::write_tsv(output[[.x]][[1]], file = paste0(output_dir, names(output)[.x], '_summary-information.tsv'))
    readr::write_tsv(output[[.x]][[2]], file = paste0(output_dir, names(output)[.x], '_results.tsv'))

    png(filename = paste0(output_dir, names(output)[.x], '-heatmaps.png'), units = 'in', width = 30, height = 10, res = 500)
    print(output[[.x]][[3]])
    dev.off()
  })
}
