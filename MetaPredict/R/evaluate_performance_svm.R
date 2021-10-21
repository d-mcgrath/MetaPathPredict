#' This function is to take in data, see what modules are present - then randomly sample the data to create incompleteness simulations and predict for all simulations, and make confusion matrices for each one. and store all results in a list.
#' @param .data Tibble object created with the read_data function - KEGG orthology annotation data for one or more bacterial genome annotation files.
#' #' @param module_vector Character vector. An optional character vector of specific KEGG Modules to scan user annotations for, and optionally to predict for. Default is NULL.
#' @importFrom magrittr "%>%"

#' @export
evaluate_svm <- function(.data, module_vector = NULL) {

  if (is.null(module_vector)) {
    module_vector = names(all_models)
  }

  responseVars <- detect_modules_svm(.data, module_vector)

  # create list of simulated increments
  ko_tibble <- purrr::map_dfr(.data, ~ {
    .x %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(genome_name) %>%
      dplyr::summarize(k_numbers = paste0(k_number, collapse = ' ')) %>%
      dplyr::as_tibble()
  }) %>%
    dplyr::mutate(k_numbers = stringr::str_split(k_numbers, pattern = ' '))

  ko_tibble <- purrr::map(seq(0.10, 1, by = 0.10), function(.prop) {
    ko_tibble %>%
      dplyr::mutate(k_numbers = purrr::map(k_numbers, ~ sample(x = .x)), # randomly shuffle the k_numbers
                    k_numbers = purrr::map(k_numbers, ~ sample(x = .x, size = .prop * length(.x))) #randomly sample the k_numbers
      ) %>%
      tidyr::unnest(cols = k_numbers) %>%
      dplyr::mutate(k_count = 1) %>%
      dplyr::group_by(genome_name, k_numbers) %>%
      dplyr::summarize(k_count = sum(k_count)) %>%
      tidyr::pivot_wider(names_from = k_numbers, values_from = k_count) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
      dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ dplyr::case_when(is.na(.x) ~ 0L,
                                                                          TRUE ~ .x))) %>%
      dplyr::bind_cols(dplyr::select(filler, -c(colnames(filler)[colnames(filler) %in% colnames(.)]))) %>%
      dplyr::select(colnames(filler)) %>%
      dplyr::relocate(colnames(filler)) %>%
      predict(caret::preProcess(., method = c('center', 'scale')), .)
  }) %>%
    purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))


  predictions <- purrr::map(ko_tibble, function(.sim_data) {
    purrr::imap_dfc(all_models[module_vector], ~
                  named_predict(.y, .x, .sim_data)
    )
  })

  confusion_matrices <- purrr::map(predictions, function(.iter) {
    purrr::map2(.iter, dplyr::select(responseVars$pres_abs_tbl, -genome_name), ~
                  if (all(colnames(.x) == colnames(.y))) {
                    caret::confusionMatrix(factor(.y, levels = c('no', 'yes')),
                                           factor(.x, levels = c('no', 'yes')),
                                           positive = 'yes')
                  } else {
                    cli::cli_alert_danger('Error: column names of predictions do not match column names of response variables.')
                    return(NULL)
                  }
    )
  })

  performance_metrics <- purrr::map(confusion_matrices, function(.iter) {
    purrr::imap_dfr(.iter, ~
                      .x %>%
                      purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
                      purrr::flatten() %>%
                      dplyr::as_tibble() %>%
                      dplyr::mutate(model_name = .y, .before = 1) %>%
                      dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1)
    )
  })

  confusion_tibble <- purrr::map(confusion_matrices, function(.iter) {
    purrr::imap_dfr(.iter, ~
                      .x$table %>%
                      broom::tidy() %>%
                      suppressWarnings() %>%
                      dplyr::mutate(model_name = .y, .before = 1)
    )
  })

  return(list(performance_metrics = performance_metrics, confusion_tibble = confusion_tibble))
}



named_predict <- function(name, model, data) {
  z <- predict(model, data)
  names(z) <- name
  return(z)
}































