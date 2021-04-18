#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param userData Created with the read_data function - KEGG Orthology data for one or more bacteria/archaea
#' @param output_dir The full or relative path to an output directory where result and summary output will be saved as TSV flatfiles
#' @param moduleVector An optional vector of specific KEGG Modules to scan user annotations and calculate probabilities for
#' @importFrom magrittr "%>%"

#' @export
MetaPredict <- function(userData, output_dir = NULL, moduleVector = NULL, strict = FALSE) {
  cli::cli_h1('Starting MetaPredict')
  cli::cli_alert_info('Formatting input data...')

  if (missing(userData)) {
    cli::cli_alert_danger('Input userData not detected in global environment. Make sure you have run read_genome_data() or read_metagenome_data() on your data and have listed your input data object for the userData argument.')
    stop()
  }

  cli::cli_alert_info('Reconstructing metabolic pathways and performing prediction calculations...')

  gene_counts.list <- get_gene_counts_from(userData)
  results <- map_modules_to(gene_counts.list, userData, strict = strict)# this function should be split into 2 parts
  results <- get_parameters(results, moduleVector = moduleVector, strict = strict)
  results <- get_posteriors(results, strict = strict)

  if (unique(userData$data_type) == 'metagenome') {
    results <- results %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(module_class, step) %>%
      dplyr::mutate(metagenome_name = unique(userData$metagenome_name)) # needs updating - to allow for multiple metagenomes as input...!

    summary_information <- summarize_metagenome_output(results)
    #heatmaps <- create_heatmaps_from(results, userData, metagenome = TRUE)
    output <- list(list('summary' = summary_information, 'full_results' = results)) %>% #, heatmaps)) %>%
      purrr::set_names(unique(userData$metagenome_name))

    if (is.null(output_dir)) {
      cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and KEGG Module probability calculations. Output is in a list.')
      cli::cli_alert_info('Enter View(results[[(x)]][[1]]) to look at summary information, View(results[[(x)]][[2]] for the full results.') #, and print(results[[(x)]][[3]]) for pathway completion heatmaps.')
      cli::cli_alert_warning('Make sure to change (x) to a number between 1 and x, with x equal to the total amount of input genome annotations.')
      cli::cli_alert_info("To save results, run the following command: save_results(output_dir = '/path/to/output/directory'). Note: the output directory will be created if it does not exist.")
      cli::cli_h1('All done.')
      return(output)
    }
  } else {
    summary_information <- summarize_genome_output(results)
    #heatmaps <- purrr::map(results, ~ {create_heatmaps_from(.x, userData)})
    output <- purrr::transpose(list('summary' = summary_information, 'full_results' = results)) #, heatmaps))

    if (is.null(output_dir)) {
      cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and reaction probability calculations. Output is in a list.')
      cli::cli_alert_info('Enter View(results[[(x)]][[1]]) to look at summary information, View(results[[(x)]][[2]] for the full results.') #, and print(results[[(x)]][[3]]) for pathway completion heatmaps.')
      cli::cli_alert_warning('Make sure to change (x) to a number from 1 to x, with x equal to the total amount of input genome/metagenome inputs')
      cli::cli_alert_info("To save results, run the following command: save_results(output_dir = '/path/to/output/directory/'). Note: the output directory will be created if it does not exist.")
      cli::cli_h1('All done.')
      return(output)
    }
  }
  if (!is.null(output_dir)) {
    save_results(output, output_dir)
    cli::cli_alert_success('Finished KEGG metabolic pathway reconstruction and reaction probability calculations. Output is in directory: {output_dir}')
    cli::cli_h1('All done.')
    return(output)
  }
}



## ADD GENOME/METAGENOME NAME AS A COLUMN
summarize_genome_output <- function(.data) {
  .data <- purrr::map(seq_along(.data), ~ {
    .data[[.x]] <- .data[[.x]] %>%
      dplyr::filter(!(duplicated(step))) %>%
      dplyr::group_by(module) %>%
      dplyr::add_tally(module_step_present == TRUE, name = 'steps_present') %>%
      dplyr::add_tally(length(module), name = 'module_length') %>%
      dplyr::add_tally(module_step_present == FALSE, name = 'predicted') %>%
      dplyr::add_tally(probability >= 0.90, name = 'p_greater_90') %>%
      dplyr::select(taxonomy, taxonomy_used, module_name, module_class, module, steps_present, module_length,
             predicted, p_greater_90) %>%
      dplyr::distinct() %>%
      dplyr::mutate('Module steps present' = paste(steps_present, module_length, sep = '/'),
             'Module steps predicted' = paste(predicted, module_length, sep = '/'),
             'Predictions (P > 0.90)' = paste(p_greater_90, module_length, sep = '/'),
             'Predicted completeness (P > 0.90)' = paste(steps_present+p_greater_90, module_length, sep = '/')) %>%
      dplyr::rename('Module name' = module_name, 'Module class' = module_class, 'Module' = module, 'Taxonomy' = taxonomy) %>%
      dplyr::select(-c(steps_present, module_length, p_greater_90, predicted))
  })
  return(.data)
}



summarize_metagenome_output <- function(.data) {
  .data %>%
    dplyr::filter(!(duplicated(step))) %>%
    dplyr::group_by(module, taxonomy_used) %>%
    dplyr::add_tally(module_step_present == TRUE, name = 'steps_present') %>%
    dplyr::add_tally(length(module), name = 'module_length') %>%
    dplyr::add_tally(module_step_present == FALSE, name = 'predicted') %>%
    dplyr::add_tally(probability >= 0.90, name = 'p_greater_90') %>%
    dplyr::select(taxonomy, taxonomy_used, module_name, module_class, module, steps_present, module_length,
           predicted, p_greater_90) %>%
    dplyr::distinct() %>%
    dplyr::mutate('Module steps present' = paste(steps_present, module_length, sep = '/'),
           'Module steps predicted' = paste(predicted, module_length, sep = '/'),
           'Predictions (P > 0.90)' = paste(p_greater_90, module_length, sep = '/'),
           'Predicted completeness (P > 0.90)' = paste(steps_present + p_greater_90, module_length, sep = '/')) %>%
    dplyr::rename('Module name' = module_name, 'Module class' = module_class, 'Module' = module, 'Taxonomy' = taxonomy) %>%
    dplyr::select(-c(steps_present, module_length, p_greater_90, predicted))
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

    #png(filename = paste0(output_dir, names(output)[.x], '-heatmaps.png'), units = 'in', width = 30, height = 10, res = 500)
    #print(output[[.x]][[3]])
    #dev.off()
  })
}
