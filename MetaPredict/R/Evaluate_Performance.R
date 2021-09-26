#' This function is to take in data, see what modules are present - then randomly sample the data to create incompleteness simulations and predict for all simulations, and make confusion matrices for each one. and store all results in a list.
#' @param .data Tibble object created with the read_data function - KEGG orthology annotation data for one or more bacterial genome annotation files.
#' #' @param moduleVector Character vector. An optional character vector of specific KEGG Modules to scan user annotations for, and optionally to predict for. Default is NULL.
#' @importFrom magrittr "%>%"

#' @export
evaluate_dataset <- function(.data, moduleVector = NULL) {

  if (is.null(moduleVector)) {
    moduleVector = names(all_models)
  }

  responseVars <- detect_modules(.data, moduleVector)

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
      as.matrix()
  }) %>%
    purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))

  #for the complete data and each set of simulated incomplete data: predict presence/absence of each KEGG module
  predictions <- purrr::map(ko_tibble, function(.sim_data) {

    purrr::map(all_models[moduleVector], ~
                 predict(.x$glmnet.fit,
                         s = .x$lambda.1se,
                         newx = .sim_data,
                         type = 'class')) %>%
      purrr::map_dfc(~ .x) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
      dplyr::rename_with(~ moduleVector)
  })

  confusion_matrices <- purrr::map(predictions, function(.iter) {
    purrr::map2(.iter, dplyr::select(responseVars$pres_abs_tbl, -genome_name), ~
                  if (all(colnames(.x) == colnames(.y))) {
                    caret::confusionMatrix(factor(.y, levels = c('0', '1')),
                                           factor(.x, levels = c('0', '1')),
                                           positive = '1')
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
                      dplyr::mutate(model_name = .y, .before = 1)
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







#' @export
evaluate_model_testdata <- function(responseVars, ko_tibble, moduleVector = NULL) {

  if (is.null(moduleVector)) {
    moduleVector = names(responseVars)
  }

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
      as.matrix()
  }) %>%
    purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))

  #for the complete data and each set of simulated incomplete data: predict presence/absence of each KEGG module
  predictions <- purrr::map(ko_tibble, function(.sim_data) {

    purrr::map(all_models[moduleVector], ~
                 predict(.x$glmnet.fit,
                         s = .x$lambda.1se,
                         newx = .sim_data,
                         type = 'class')) %>%
      purrr::map_dfc(~ .x) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
      dplyr::rename_with(~ moduleVector)
  })

  confusion_matrices <- purrr::map(predictions, function(.iter) {
    purrr::map2(.iter, responseVars, ~
                  if (all(colnames(.x) == colnames(.y))) {
                    caret::confusionMatrix(factor(.y, levels = c('0', '1')),
                                           factor(.x, levels = c('0', '1')),
                                           positive = '1')
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
                      dplyr::mutate(model_name = .y, .before = 1)
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



#' @export
evaluate_model_testdata_dt <- function(responseVars, ko_tibble, moduleVector = NULL) {

  if (is.null(moduleVector)) {
    moduleVector <- names(responseVars)
  }

  ko_tibble <- purrr::map(seq(0.10, 1, by = 0.10), function(.prop) {
    ko_tibble %>%
      dplyr::mutate(k_numbers = purrr::map(k_numbers, ~ sample(x = .x)), # randomly shuffle the k_numbers
                    k_numbers = purrr::map(k_numbers, ~ sample(x = .x, size = .prop * length(.x))) #randomly sample the k_numbers
      ) %>%
      tidyr::unnest(cols = k_numbers) %>%
      dplyr::mutate(k_count = 1) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(genome_name, k_numbers) %>%
      dplyr::summarize(k_count = sum(k_count)) %>%
      tidyr::pivot_wider(names_from = k_numbers, values_from = k_count) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
      dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ dplyr::case_when(is.na(.x) ~ 0L,
                                                                          TRUE ~ .x))) %>%
      dplyr::as_tibble() %>%
      dplyr::bind_cols(dplyr::select(filler, -c(colnames(filler)[colnames(filler) %in% colnames(.)]))) %>%
      dplyr::select(colnames(filler)) %>%
      dplyr::relocate(colnames(filler)) %>%
      as.matrix()
  }) %>%
    purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))

  #for the complete data and each set of simulated incomplete data: predict presence/absence of each KEGG module
  predictions <- purrr::map(ko_tibble, function(.sim_data) {

    purrr::map(all_models[moduleVector], ~
                 predict(.x$glmnet.fit,
                         s = .x$lambda.1se,
                         newx = .sim_data,
                         type = 'class')) %>%
      purrr::map_dfc(~ .x) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
      dplyr::rename_with(~ moduleVector)
  })

  confusion_matrices <- purrr::map(predictions, function(.iter) {
    purrr::map2(.iter, responseVars, ~
                  if (all(colnames(.x) == colnames(.y))) {
                    caret::confusionMatrix(factor(.y, levels = c('0', '1')),
                                           factor(.x, levels = c('0', '1')),
                                           positive = '1')
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
                      dplyr::mutate(model_name = .y, .before = 1)
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



#' @export
evaluate_model_testdata_parallel <- function(responseVars, ko_tibble, moduleVector = NULL) {

  future::plan(multicore, workers = 69)
  options(future.globals.maxSize = 32505856000) #31000 * 1024 ^2


  if (is.null(moduleVector)) {
    moduleVector <- names(responseVars)
  }

  ko_tibble <- furrr::future_map(furrr::furrr_options(seed = TRUE), seq(0.10, 1, by = 0.10), function(.prop) {
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
      as.matrix()
  }) %>%
    purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))

  #for the complete data and each set of simulated incomplete data: predict presence/absence of each KEGG module
  predictions <- purrr::map(ko_tibble, function(.sim_data) {

    furrr::future_map(furrr::furrr_options(seed = TRUE), all_models[moduleVector], ~
                 predict(.x$glmnet.fit,
                         s = .x$lambda.1se,
                         newx = .sim_data,
                         type = 'class')) %>%
      purrr::map_dfc(~ .x) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.double(.x))) %>%
      dplyr::rename_with(~ moduleVector)
  })

  confusion_matrices <- furrr::future_map(furrr::furrr_options(seed = TRUE), predictions, function(.iter) {
    purrr::map2(.iter, responseVars, ~
                  if (all(colnames(.x) == colnames(.y))) {
                    caret::confusionMatrix(factor(.y, levels = c('0', '1')),
                                           factor(.x, levels = c('0', '1')),
                                           positive = '1')
                  } else {
                    cli::cli_alert_danger('Error: column names of predictions do not match column names of response variables.')
                    return(NULL)
                  }
    )
  })

  performance_metrics <- furrr::future_map(furrr::furrr_options(seed = TRUE), confusion_matrices, function(.iter) {
    purrr::imap_dfr(.iter, ~
                      .x %>%
                      purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
                      purrr::flatten() %>%
                      dplyr::as_tibble() %>%
                      dplyr::mutate(model_name = .y, .before = 1)
    )
  })

  confusion_tibble <- furrr::future_map(furrr::furrr_options(seed = TRUE), confusion_matrices, function(.iter) {
    purrr::imap_dfr(.iter, ~
                      .x$table %>%
                      broom::tidy() %>%
                      dplyr::mutate(model_name = .y, .before = 1)
    )
  })

  return(list(performance_metrics = performance_metrics, confusion_tibble = confusion_tibble))
}


#' @export
test_fn <- function(responseVars, ko_tibble, moduleVector = NULL, filler = filler) {

  future::plan(future::multicore, workers = 69)
  options(future.globals.maxSize = 32505856000) #31000 * 1024 ^2

  if (is.null(moduleVector)) {
    moduleVector <- names(responseVars)
  }

  ko_tibble <- furrr::future_map(furrr::furrr_options(seed = TRUE), seq(0.10, 1, by = 0.10), function(.prop) {
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
      as.matrix()
  }) #%>%
  #  purrr::set_names(nm = paste0('prop.', seq(10, 100, by = 10)))

  return(ko_tibble)
}



