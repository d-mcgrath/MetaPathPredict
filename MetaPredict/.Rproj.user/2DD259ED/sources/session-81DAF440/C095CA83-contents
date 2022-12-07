
# Various helper functions ------------------------------------------------

reconstruct_modules <- function(from, module_vector = module_vector, module_detect_type = module_detect_type) {
  if (module_detect_type == 'extract') {
    reconstructed <- extract_modules(from, .modules = module_vector)
  } else {
    reconstructed <- detect_modules(from, .modules = module_vector)
  }
  return(reconstructed)
}



check_data <- function(data, module_vector) {
  if (missing(data)) {
    cli::cli_alert_danger('Input data not detected in global environment. Make sure you have run read_data() on your data and have listed your input data object for the data argument.')
    stop()
  }

  if (any(!(module_vector %in% all_model_names))) {
    bad_modules <- module_vector[!(module_vector %in% all_model_names)]
    cli::cli_alert_danger('The following module(s) are improperly named or not available at this time: {bad_modules}. Use available_modules() to check which KEGG Modules are currently available for predictions.')
    stop()
  }
  return(data)
}



#' @export
named_predict <- function(model_name, data, database_reference) {

  model <- database_reference |>
    dplyr::filter(model_name == model_name) |>
    dplyr::collect() |>
    dplyr::pull(raw_model) |>
    unserialize_model()

  prediction <- predict(model, data)

  names(prediction) <- model_name
  return(prediction)
}



# function to properly index and then unserialize a raw model "blob"
unserialize_model = function(.data) {
  unserialize(.data[[1]])
}



#' @export
create_kegg_matrix <- function(.data) {
  feature_table <- purrr::map(.data, ~ {
    .x |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(copy_number = 1) |>
      dplyr::group_by(k_number) |>
      dplyr::summarize(copy_number = sum(copy_number)) |>
      tidyr::pivot_wider(names_from = k_number, values_from = copy_number) |>
      dplyr::as_tibble()
  }) |>
    purrr::map_dfr(~ .x) |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(.x))) |>
    dplyr::as_tibble()

  feature_table <- feature_table |>
    dplyr::bind_cols(
      dplyr::select(filler,
                    -c(colnames(filler)[colnames(filler) %in% colnames(feature_table)]))) |>
    na_to_zero() |>
    dplyr::select(colnames(filler)) |>
    dplyr::relocate(colnames(filler))
}



#' @export
detect_modules <- function(.data, .modules) {
  kegg_rules <- kegg_rules |> dplyr::filter(module %in% .modules)

  detected_modules <- kegg_rules |>
    dplyr::bind_cols(
      purrr::imap_dfc(
        .data |> purrr::set_names(purrr::map_chr(.data, ~ unique(.x |> pull(genome_name)))), ~
          kegg_rules |>
          dplyr::mutate(!!.y := purrr::map_lgl(rule, function(.rule) {
            any(.rule %in% .x$k_number)})) |>
          dplyr::select(!!.y))) |>
    dplyr::group_by(module, step) |>
    dplyr::summarize(dplyr::across(tidyselect::where(is.logical), ~ all(.x))) |>
    dplyr::summarize(dplyr::across(tidyselect::where(is.logical), ~ all(.x))) |>
    dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) |>
    tidyr::pivot_longer(
      cols = 2:last_col(),
      names_to = 'genome_name',
      values_to = 'prediction') |>
    tidyr::pivot_wider(
      names_from = module,
      values_from = prediction)

  return(list(pres_abs_tbl = detected_modules))
}





# depr_detect_modules <- function(.data, .modules) {
#   patt.kegg_modules <- dplyr::filter(patt.kegg_modules, module %in% .modules)
#   all_kegg_modules <- dplyr::filter(all_kegg_modules, module %in% .modules)
#
#   detected.modules <- purrr::map(.data, ~ {
#     .x %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::group_by(genome_name) %>%
#       dplyr::summarize(k_number = paste0(k_number, collapse = ' ')) %>%
#       dplyr::as_tibble()
#   }) %>%
#     purrr::map_dfr(~ .x) %>%
#     dplyr::bind_cols(
#       purrr::set_names(.$k_number, nm = .$genome_name) %>%
#         purrr::map_dfc(~ stringr::str_detect(string = .x, pattern = patt.kegg_modules$k_numbers)) %>%
#         dplyr::mutate(module_step = patt.kegg_modules$step, .before = 1) %>%
#         dtplyr::lazy_dt() %>%
#         dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
#         dplyr::group_by(module_step) %>%
#         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>% # aggregate sums of multi-step steps
#         dplyr::mutate(module_name = stringr::str_replace(module_step, '(M\\d{5}).*', '\\1'),
#                       n = all_kegg_modules[['n']]) %>%
#         dplyr::relocate(c(module_name, n), .before = module_step) %>%
#         dplyr::select(-module_step) %>%
#         dplyr::group_by(module_name) %>%
#         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
#         dplyr::mutate(across(3:dplyr::last_col(), ~ dplyr::case_when(.x == n ~ 1L,
#                                                                      .x < n ~ 0L))) %>% # this does not account for multiple ways to perform one step; it requires all alternate step possibilities to be present to count the pathway as complete. need to revise this
#         dplyr::select(-n) %>%
#         dplyr::as_tibble() %>%
#         tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = 'temp', values_to = 'temp_values') %>%
#         tidyr::pivot_wider(names_from = 'module_name', values_from = 'temp_values')
#     ) %>%
#     dplyr::select(-c(temp, k_number))
#   return(list(pres_abs_tbl = detected.modules))
# }




summarize_results <- function(.recon, .pred, .module_metadata) {
  result <- .pred |>
    tidyr::pivot_longer(2:dplyr::last_col(), names_to = 'module', values_to = 'values') |>
    tidyr::pivot_wider(names_from = genome_name, values_from = values) |>
    dplyr::left_join(.module_metadata, by = 'module') |>
    dplyr::relocate(c(module, module_name, module_class, n), .before = 1)

  return(result)
}




summarize_reconstruction <- function(.recon, .module_metadata) {
  result <- .recon |>
    tidyr::pivot_longer(2:dplyr::last_col(), names_to = 'module', values_to = 'values') |>
    tidyr::pivot_wider(names_from = genome_name, values_from = values) |>
    dplyr::left_join(.module_metadata, by = 'module') |>
    dplyr::relocate(c(module, module_name, module_class, n), .before = 1)
  return(result)
}




save_results <- function(.results, output_dir, output_prefix = NULL, overwrite = FALSE) {

  .results = list(.results@summary, .results@module_reconstructions,
                  .results@module_predictions, .results@module_extractions)

  if (!dir.exists(output_dir)) {
    cli::cli_alert_warning('Creating output directory {output_dir}')
    dir.create(output_dir)
  }
  if (!stringr::str_detect(output_dir, '.*/$')) {
    output_dir <- sub('(.*)', '\\1\\/', output_dir, perl = TRUE)
  }

  if (!is.null(output_prefix)) {
    output_prefix = paste0(output_prefix, '_')
  }

  if (length(.results) == 3) {

    if (all(!(paste0(output_prefix, c('summary', 'module_reconstructions', 'module_predictions'), '.tsv') %in%
              list.files(path = output_dir)))) {
      purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions'), ~ {
        vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))})
    } else {
      if (overwrite == TRUE) {
        cli::cli_alert_warning('Overwriting existing files with output files.')
        purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions'), ~ {
          vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
        })
      } else {
        for (.attempt in 1:3) {
          ans <- readline(prompt = paste0('One or more files exist with the same name as the output to be saved. Overwrite these existing files? (attempt ', .attempt, ' of 3) [y/n] '))
          if (ans == 'y') {
            cli::cli_alert_warning('Overwriting existing files with output files.')
            purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions'), ~ {
              vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
            })
            break
          } else if (ans == 'n') {
            cli::cli_alert_info('Not saving any output files. To save files, run save_results() and try changing the "output_prefix" argument to a unique identifier.')
            break
          } else {
            if (.attempt < 3) {
              cli::cli_alert_danger('Please enter either "y" or "n" and press enter.')
            } else {
              cli::cli_alert_warning('Did not recognize what was entered. To save files, run save_results() and try changing the "output_prefix" argument to a unique identifier.')
            }
          }
        }
      }
    }
  } else {
    if (all(!(paste0(output_prefix, c('summary', 'module_reconstructions', 'module_predictions', 'module_extractions'), '.tsv') %in% list.files(path = output_dir)))) {
      purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions', 'module_extractions'), ~ {
        vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
      })
    } else {
      if (overwrite == TRUE) {
        cli::cli_alert_warning('Overwriting existing files with output files.')
        purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions', 'module_extractions'), ~ {
          vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
        })
      } else {
        for (.attempt in 1:3) {
          ans <- readline(prompt = paste0('One or more files exist with the same name as the output to be saved. Overwrite these existing files? (attempt ', .attempt, ' of 3) [y/n] '))
          if (ans == 'y') {
            cli::cli_alert_warning('Overwriting existing files with output files.')
            purrr::walk2(1:length(.results), c('summary', 'module_reconstructions', 'module_predictions', 'module_extractions'), ~ {
              vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
            })
            break
          } else if (ans == 'n') {
            cli::cli_alert_info('Not saving any output files. To save files, run save_results() and try changing the "output_prefix" argument to a unique identifier.')
            break
          } else {
            if (.attempt < 3) {
              cli::cli_alert_danger('Please enter either "y" or "n" and press enter.')
            } else {
              cli::cli_alert_warning('Did not recognize what was entered. To save files, run save_results() and try changing the "output_prefix" argument to a unique identifier.')
            }
          }
        }
      }
    }
  }

}




save_recon <- function(.results, output_dir, output_prefix = NULL, overwrite = FALSE) {

  .results = list(.results@summary, .results@module_reconstructions)

  if (!dir.exists(output_dir)) {
    cli::cli_alert_warning('Creating output directory {output_dir}')
    dir.create(output_dir)
  }
  if (!stringr::str_detect(output_dir, '.*/$')) {
    output_dir <- sub('(.*)', '\\1\\/', output_dir, perl = TRUE)
  }

  if (!is.null(output_prefix)) {
    output_prefix = paste0(output_prefix, '_')
  }

  if (all(!(paste0(output_prefix, c('summary', 'module_reconstructions'), '.tsv') %in% list.files(path = output_dir)))) {
    purrr::map2(1:length(.results), c('summary', 'module_reconstructions'), ~ {
      vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
    })
  } else {
    if (overwrite == TRUE) {
      cli::cli_alert_warning('Overwriting existing files with output files.')
      purrr::map2(1:length(.results), c('summary', 'module_reconstructions'), ~ {
        vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
      })
    } else {
      for (.attempt in 1:3) {
        ans <- readline(prompt = paste0('One or more files exist with the same name as the output to be saved. Overwrite these existing files? (attempt ', .attempt, ' of 3) [y/n] '))
        if (ans == 'y') {
          cli::cli_alert_warning('Overwriting existing files with output files.')
          purrr::map2(1:length(.results), c('summary', 'module_reconstructions'), ~ {
            vroom::vroom_write(.results[[.x]], file = paste0(output_dir, output_prefix, .y, '.tsv'))
          })
          break
        } else if (ans == 'n') {
          cli::cli_alert_info('Not saving any output files. To save files, run save_recon() and try changing the "output_prefix" argument to a unique identifier.')
          break
        } else {
          if (.attempt < 3) {
            cli::cli_alert_danger('Please enter either "y" or "n" and press enter.')
          } else {
            cli::cli_alert_warning('Did not recognize what was entered. To save files, run save_recon() and try changing the "output_prefix" argument to a unique identifier.')
          }
        }
      }
    }
  }
}




extract_modules <- function(.data, .modules) {
  patt.kegg_modules <- dplyr::filter(patt.kegg_modules, module %in% .modules)
  all_kegg_modules <- dplyr::filter(all_kegg_modules, module %in% .modules)

  .data_cat <- purrr::map_dfr(.data, ~ {
    .x |>
      dtplyr::lazy_dt() |>
      dplyr::group_by(genome_name) |>
      dplyr::summarize(k_number = paste0(k_number, collapse = ' ')) |>
      dplyr::as_tibble()
  })

  extractions <- purrr::map2_dfc(.data_cat$k_number, .data_cat$genome_name, ~
                                   tibble::tibble(!!.y := stringi::stri_extract_all_regex(
                                     str = .x, pattern = patt.kegg_modules$k_numbers))) |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(module = patt.kegg_modules$alt_steps,
                  group_col = dplyr::row_number()) |>
    dplyr::relocate(c(module, group_col), .before = 1) |>
    dplyr::group_by(module, group_col) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ paste0(unlist(.x[!is.na(.x)]), collapse = '_'))) |>  # for each step (or 1 step of multi-step step): glue together duplicate genes or diff genes that do the same step (e.g., in case when requirement is K01626|K03856, if both are present and 2 copies of K01626 present, it will become: K01626_K01626_K03856)
    dplyr::ungroup() |>
    dplyr::select(-group_col) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringi::stri_replace_first_regex(.x, '^$', NA_character_))) |>
    dplyr::group_by(module) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), ~ if (any(is.na(.x))) {NA_character_} else {paste0(unlist(.x[!is.na(.x)]), collapse = '_')} )) |>  # This aggregates all steps in multi-step steps that have an _A and _B alternate way of completing - and condenses each _A and _B method into a single row, with the genes separated by '_'; if any piece of an alt step is missing, it returns NA instead
    dplyr::ungroup() |>
    dplyr::mutate(module = stringr::str_replace(module, '(M\\d{5}_\\d{2}_)[:ALPHA:].*', '\\1alt')) |>
    dplyr::group_by(module) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), ~ paste0(unlist(.x[!is.na(.x)]), collapse = '|'))) |>  #combine alternate versions of the same module step with `|`
    dplyr::ungroup() |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringi::stri_replace_first_regex(.x, '^$', NA_character_))) |>
    dplyr::mutate(module = stringr::str_replace(module, '(M\\d{5}_\\d{2}).*', '\\1')) |>
    dplyr::mutate(n_steps = 1) |>
    dplyr::relocate(n_steps, .after = 1) |>
    dplyr::group_by(module) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), ~ if (any(is.na(.x))) {NA_character_} else {paste0(unlist(.x[!is.na(.x)]), collapse = '+')} )) |>  #combine multi-step module steps together with `+`; returns NA if not all steps of multi-step steps are present, or if the step of a single step step is missing
    dplyr::ungroup() |>
    dplyr::mutate(module = stringr::str_replace(module, '(M\\d{5}).*', '\\1')) |>
    dplyr::group_by(module) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), ~ paste0(.x[!is.na(.x)], collapse = ' '))) |>
    dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ stringi::stri_replace_all_regex(.x, '^$', NA_character_))) |>
    dplyr::mutate(n_steps = (stringi::stri_count_fixed(n_steps, ' ') + 1)) |>
    dplyr::as_tibble()

  module_pres_abs_int_tbl <- extractions |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(dplyr::across(3:dplyr::last_col(), ~ (stringi::stri_count_fixed(.x, ' ') + 1))) |>
    dplyr::mutate(dplyr::across(3:dplyr::last_col(), ~ dplyr::case_when(is.na(.x) ~ 0, TRUE ~ .x))) |>
    dplyr::mutate(dplyr::across(3:dplyr::last_col(), ~ dplyr::case_when(.x == n_steps ~ 1L, TRUE ~ 0L))) |>
    dplyr::select(-n_steps) |>
    dplyr::as_tibble() |>
    tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = 'genome_name', values_to = 'temp_values') |>
    tidyr::pivot_wider(names_from = 'module', values_from = 'temp_values')

  return(list(ko_extractions = extractions, pres_abs_tbl = module_pres_abs_int_tbl))
}




# detect_modules.corrected <- function(.data, .modules) {
#   patt.kegg_modules <- dplyr::filter(patt.kegg_modules, module %in% .modules)
#   all_kegg_modules <- dplyr::filter(all_kegg_modules, module %in% .modules)
#   module_metadata <- dplyr::filter(module_metadata, module %in% .modules)
#
#   detected.modules <- purrr::map_dfr(.data, ~ {
#     .x %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::group_by(genome_name) %>%
#       dplyr::summarize(k_number = paste0(k_number, collapse = ' ')) %>%
#       dplyr::as_tibble()
#   }) %>%
#     dplyr::bind_cols(
#       purrr::set_names(.$k_number, nm = .$genome_name) %>%
#         purrr::map_dfc(~ stringi::stri_detect_regex(str = .x, pattern = patt.kegg_modules$k_numbers)) %>%
#         dplyr::mutate(module_step = patt.kegg_modules$step, .before = 1) %>%
#         dtplyr::lazy_dt() %>%
#         dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
#         dplyr::group_by(module_step) %>%
#         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>% # aggregate sums of multi-step steps
#         dplyr::mutate(n = all_kegg_modules[['n']]) %>%
#         dplyr::relocate(n, .after = 1) %>%
#         dplyr::mutate(dplyr::across(3:dplyr::last_col(), ~ dplyr::if_else(n == .x, 1, 0))) %>%
#         dplyr::select(-n) %>%
#         dplyr::mutate(module_step = stringr::str_replace(module_step, '(M\\d{5}_\\d{2}_)[:ALPHA:]', '\\1alt')) %>%
#         dplyr::group_by(module_step) %>%
#         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ dplyr::if_else(.x >= 1, 1, 0))) %>%
#
#         dplyr::mutate(module_name = stringr::str_replace(module_step, '(M\\d{5}).*', '\\1')) %>%
#         dplyr::relocate(module_name, .before = 1) %>%
#         dplyr::select(-module_step) %>%
#         dplyr::group_by(module_name) %>%
#         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
#         dplyr::mutate(n = module_metadata[['n']]) %>%
#         dplyr::relocate(n, .after = 1) %>%
#         dplyr::mutate(dplyr::across(3:dplyr::last_col(), ~ dplyr::case_when(.x == n ~ 1L,
#                                                                             .x < n ~ 0L))) %>% # this does not account for multiple ways to perform one step; it requires all alternate step possibilities to be present to count the pathway as complete. need to revise this
#         dplyr::select(-n) %>%
#         dplyr::as_tibble() %>%
#         tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = 'temp', values_to = 'temp_values') %>%
#         tidyr::pivot_wider(names_from = 'module_name', values_from = 'temp_values')
#     ) %>%
#     dplyr::select(-c(temp, k_number))
#   return(list(pres_abs_tbl = detected.modules))
# }



#' detect_modules_svm <- function(.data, .modules) {
#'   patt.kegg_modules <- dplyr::filter(patt.kegg_modules, module %in% .modules)
#'   all_kegg_modules <- dplyr::filter(all_kegg_modules, module %in% .modules)
#'
#'   detected.modules <- purrr::map(.data, ~ {
#'     .x %>%
#'       dtplyr::lazy_dt() %>%
#'       dplyr::group_by(genome_name) %>%
#'       dplyr::summarize(k_number = paste0(k_number, collapse = ' ')) %>%
#'       dplyr::as_tibble()
#'   }) %>%
#'     purrr::map_dfr(~ .x) %>%
#'     dplyr::bind_cols(
#'       purrr::set_names(.$k_number, nm = .$genome_name) %>%
#'         purrr::map_dfc(~ stringr::str_detect(string = .x, pattern = patt.kegg_modules$k_numbers)) %>%
#'         dplyr::mutate(module_step = patt.kegg_modules$step, .before = 1) %>%
#'         dtplyr::lazy_dt() %>%
#'         dplyr::mutate(dplyr::across(2:dplyr::last_col(), ~ as.integer(.x))) %>%
#'         dplyr::group_by(module_step) %>%
#'         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>% # aggregate sums of multi-step steps
#'         dplyr::mutate(module_name = stringr::str_replace(module_step, '(M\\d{5}).*', '\\1'),
#'                       n = all_kegg_modules[['n']]) %>%
#'         dplyr::relocate(c(module_name, n), .before = module_step) %>%
#'         dplyr::select(-module_step) %>%
#'         dplyr::group_by(module_name) %>%
#'         dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) %>%
#'         dplyr::mutate(across(3:dplyr::last_col(), ~ dplyr::case_when(.x == n ~ 'yes',
#'                                                                      .x < n ~ 'no'))) %>% # this does not account for multiple ways to perform one step; it requires all alternate step possibilities to be present to count the pathway as complete. need to revise this
#'         dplyr::select(-n) %>%
#'         dplyr::as_tibble() %>%
#'         tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = 'temp', values_to = 'temp_values') %>%
#'         tidyr::pivot_wider(names_from = 'module_name', values_from = 'temp_values')
#'     ) %>%
#'     dplyr::select(-c(temp, k_number))
#'   return(list(pres_abs_tbl = detected.modules))
#' }
