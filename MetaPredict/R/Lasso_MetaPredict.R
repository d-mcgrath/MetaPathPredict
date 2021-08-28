#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param .data Created with the read_data function - KEGG Orthology data for one or more bacteria/archaea
#' @param output_dir The full or relative path to an output directory where result and summary output will be saved as TSV flatfiles
#' @param moduleVector An optional vector of specific KEGG Modules to scan user annotations and calculate probabilities for
#' @importFrom magrittr "%>%"

#' @export
MetaPredict <- function(.data, output_dir = NULL, moduleVector = NULL, predict_models = TRUE, predict_type = c('response', 'class')) {
  cli::cli_h1('Starting MetaPredict')
  #cli::cli_alert_info('Formatting input data...')

  if (missing(.data)) {
    cli::cli_alert_danger('Input data not detected in global environment. Make sure you have run read_data() on your data and have listed your input data object for the .data argument.')
    stop()
  }

  cli::cli_alert_info('Reconstructing KEGG modules presence/absence...')

  if (is.null(moduleVector)) {
    reconstructed <- detect_modules(.data, .modules = names(all_models))

  if (predict_models == FALSE) {
    summary <- summarize_recon(.recon = reconstructed, .module_metadata = module_metadata)
    results <- list(summary = summary, module_reconstructions = reconstructed)
    if (!is.null(output_dir)) {
      cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
      save_recon(results, output_dir)
    } else {
      cli::cli_alert_success('All done.')
      cli::cli_alert_info('To save results, use save_recon().')
    }
    return(results)
  }
  cli::cli_alert_info('Performing prediction calculations...')

  .data <- create_kegg_matrix(.data) %>%
    as.matrix()

  predictions <- purrr::map(all_models, ~ predict(.x$glmnet.fit,
                                 s = .x$lambda.1se,
                                 newx = .data,
                                 type = predict_type)) %>%
    purrr::map_dfc(~ .x) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
    dplyr::rename_with(~ names(all_models)) %>%
    dplyr::mutate(genome_name = reconstructed$genome_name, .before = 1)

  #return the results
  for (.column in colnames(reconstructed)[2:length(reconstructed)]) {
    predictions[reconstructed[, .column] == 1L, .column] <- NA_real_
  }

  summary <- summarize_results(.recon = reconstructed, .pred = predictions, .module_metadata = module_metadata)
  results <- list(summary = summary, module_reconstructions = reconstructed, module_predictions = predictions)

  if (!is.null(output_dir)) {
    cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
    save_results(results, output_dir)
  } else {
    cli::cli_alert_success('All done.')
    cli::cli_alert_info('To save results, use save_results().')
  }

  return(results)


  } else {

    reconstructed <- detect_modules(.data, .modules = moduleVector)

    if (predict_models == FALSE) {
      summary <- summarize_recon(.recon = reconstructed, .module_metadata = module_metadata)
      results <- list(summary = summary, module_reconstructions = reconstructed)
      if (!is.null(output_dir)) {
        cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
        save_recon(results, output_dir)
      } else {
        cli::cli_alert_success('All done.')
        cli::cli_alert_info('To save results, use save_recon().')
      }
      return(results)
    }

    cli::cli_alert_info('Performing prediction calculations...')

    .data <- create_kegg_matrix(.data) %>%
      as.matrix()

    predictions <- purrr::map(all_models[moduleVector], ~ predict(.x$glmnet.fit,
                                                                  s = .x$lambda.1se,
                                                                  newx = .data,
                                                                  type = predict_type)) %>%
      purrr::map_dfc(~ .x) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
      dplyr::rename_with(~ names(all_models[moduleVector])) %>%
      dplyr::mutate(genome_name = reconstructed$genome_name, .before = 1)


    #return the results
    for (.column in colnames(reconstructed)[2:length(reconstructed)]) {
      predictions[reconstructed[, .column] == 1L, .column] <- NA_real_
    }

    summary <- summarize_results(.recon = reconstructed, .pred = predictions, .module_metadata = module_metadata)
    results <- list(summary = summary, module_reconstructions = reconstructed, module_predictions = predictions)

    if (!is.null(output_dir)) {
      cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
      save_results(results, output_dir)
    } else {
      cli::cli_alert_success('All done.')
      cli::cli_alert_info('To save results, use save_results().')
    }


    return(results)
  }
}


#' @export
create_kegg_matrix <- function(.data) {
  purrr::map(.data, ~ {
    .x %>%
      dplyr::mutate(copy_number = 1) %>%
      dplyr::group_by(k_number) %>%
      dplyr::summarize(copy_number = sum(copy_number), .groups = 'drop') %>%
      tidyr::pivot_wider(names_from = k_number, values_from = copy_number)
  }) %>%
    purrr::map_dfr(~ .x) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(.x))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_when(is.na(.x) ~ 0L,
                                                                TRUE ~ .x))) %>%
    dplyr::bind_cols(dplyr::select(filler, -c(colnames(filler)[colnames(filler) %in% colnames(.)]))) %>%
    dplyr::select(colnames(filler)) %>%
    dplyr::relocate(colnames(filler))
}


#' @export
detect_modules <- function(.data, .modules) {
  patt.kegg_modules <- dplyr::filter(patt.kegg_modules, module %in% .modules)
  all_kegg_modules <- dplyr::filter(all_kegg_modules, module %in% .modules)

  detected.modules <- purrr::map(.data, ~ {
    .x %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(genome_name) %>%
      dplyr::summarize(k_number = paste0(k_number, collapse = ' ')) %>%
      dplyr::as_tibble()
  }) %>%
    purrr::map_dfr(~ .x) %>%
    dplyr::bind_cols(
      purrr::set_names(.$k_number, .$genome_name) %>%
      purrr::map_dfc(~ stringr::str_detect(string = .x, pattern = patt.kegg_modules$k_numbers)) %>%
        dplyr::mutate(module_step = patt.kegg_modules$step, .before = 1) %>%
        dtplyr::lazy_dt() %>%
        dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
        dplyr::group_by(module_step) %>%
        dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
        dplyr::mutate(module_name = stringr::str_replace(module_step, '(M\\d{5}).*', '\\1'),
                      n = all_kegg_modules[['n']]) %>%
        dplyr::relocate(c(module_name, n), .before = module_step) %>%
        dplyr::select(-module_step) %>%
        dplyr::group_by(module_name) %>%
        dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
        dplyr::mutate(across(3:dplyr::last_col(), ~ dplyr::case_when(.x == n ~ 1L,
                                                .x < n ~ 0L))) %>%
        dplyr::select(-n) %>%
        dplyr::as_tibble() %>%
        tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = 'temp', values_to = 'temp_values') %>%
        tidyr::pivot_wider(names_from = 'module_name', values_from = 'temp_values')
    ) %>%
    dplyr::select(-c(temp, k_number))
  return(detected.modules)
}



#' @export
summarize_results <- function(.recon, .pred, .module_metadata) {
  result <- .recon %>%
    dplyr::select(-genome_name) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.x))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_when(.x == '1' ~ 'Present',
                                                                        .x == '0' ~ NA_character_)))

  .pred <- dplyr::mutate(.pred, dplyr::across(dplyr::everything(), ~ as.character(.x)))

  for (.column in colnames(result))  {
    result[is.na(result[, .column]), .column] <- .pred[is.na(result[, .column]), .column]
    }

  result <- result %>%
    dplyr::mutate(genome_name = .recon$genome_name, .before = 1) %>%
    tidyr::pivot_longer(2:dplyr::last_col(), names_to = 'module', values_to = 'values') %>%
    tidyr::pivot_wider(names_from = genome_name, values_from = values) %>%
    dplyr::left_join(.module_metadata, by = 'module') %>%
    dplyr::relocate(c(module, module_name, module_class), .before = 1)

  return(result)
}



#' @export
summarize_recon <- function(.recon, .module_metadata) {
  result <- .recon %>%
    dplyr::select(-genome_name) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.x))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_when(.x == '1' ~ 'Present',
                                                                        .x == '0' ~ 'Absent')))

  result <- result %>%
    dplyr::mutate(genome_name = .recon$genome_name, .before = 1) %>%
    tidyr::pivot_longer(2:dplyr::last_col(), names_to = 'module', values_to = 'values') %>%
    tidyr::pivot_wider(names_from = genome_name, values_from = values) %>%
    dplyr::left_join(.module_metadata, by = 'module') %>%
    dplyr::relocate(c(module, module_name, module_class), .before = 1)

  return(result)
}



#' @export
save_results <- function(.results, output_dir) {
  if (!dir.exists(output_dir)) {
    cli::cli_alert_warning('Creating output directory {output_dir}')
    dir.create(output_dir)
  }
  if (!stringr::str_detect(output_dir, '.*/$')) {
    output_dir <- sub('(.*)', '\\1\\/', output_dir, perl = TRUE)
  }
  purrr::map2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions'), ~ {
    readr::write_tsv(.results[[.x]], file = paste0(output_dir, .y, '.tsv'))
  })
}



#' @export
save_recon <- function(.results, output_dir) {
  if (!dir.exists(output_dir)) {
    cli::cli_alert_warning('Creating output directory {output_dir}')
    dir.create(output_dir)
  }
  if (!stringr::str_detect(output_dir, '.*/$')) {
    output_dir <- sub('(.*)', '\\1\\/', output_dir, perl = TRUE)
  }
  purrr::map2(1:length(.results), c('summary', 'module_reconstructions'), ~ {
    readr::write_tsv(.results[[.x]], file = paste0(output_dir, .y, '.tsv'))
  })
}
