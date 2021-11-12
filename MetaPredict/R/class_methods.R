#' @importFrom magrittr "%>%"

# initialize methods ------------------------------------------------------

setMethod(f = 'initialize',
          signature = 'MetaReconstructInput',
          definition = function(.Object, ...) {
            .Object <- callNextMethod()
            validObject(.Object)
            .Object
          })



setMethod(f = 'initialize',
          signature = 'MetaPredictInput',
          definition = function(.Object, ...) {
            .Object <- callNextMethod()
            validObject(.Object)
            .Object
          })



setMethod(f = 'initialize',
          signature = 'MetaReconstructResult',
          definition = function(.Object, ...) {
            .Object <- callNextMethod()
            validObject(.Object)
            .Object
          })



setMethod(f = 'initialize',
          signature = 'MetaPredictResult',
          definition = function(.Object, ...) {
            .Object <- callNextMethod()
            validObject(.Object)
            .Object
          })





# generic function definitions --------------------------------------------

setGeneric('metapredict', function(object, ...) {
  standardGeneric('metapredict')
})



#setGeneric('save_results', function(object, ...) {
#  standardGeneric('save_results')
#})






# generic function methods ------------------------------------------------

#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param object Tibble object created with the read_data function - KEGG orthology annotation data for one or more bacterial genome annotation files.
#' @param module_vector Character vector. An optional character vector of specific KEGG Modules to scan user annotations for, and optionally to predict for. Default is a character vector containing identifiers for all available modules.
#' @param module_detect_type Character. Method 'extract' will give detailed output component containing all extracted KO gene identifiers associated with KEGG modules. Method 'detect' will only detect presence/absence of modules and will not provide detailed output in results.
#' @param output_dir Character. The full or relative path to an output directory where result and summary output will be saved as TSV flatfiles
#' @param output_prefix Character. Optional string to prefix to MetaPredict output files. Default is NULL, resulting in output files with default names.
#' @param overwrite Logical. If TRUE, output files from running the MetaPredict function that have the same name as existing files in the output directory will overwrite those existing files. Default is FALSE.

setMethod(f = 'metapredict',
          signature = 'MetaReconstructInput',
          definition = function(object, module_vector = names(all_models), #predict_models = TRUE,
                                module_detect_type = c('extract', 'detect'),
                                output_dir = NULL, output_prefix = NULL, overwrite = FALSE) {

            stopifnot(is(object, 'MetaReconstructInput'))

            cli::cli_h1('Starting MetaPredict')
            module_detect_type <- match.arg(module_detect_type)
            object <- check_data(object, module_vector = module_vector)

            cli::cli_alert_info('Reconstructing KEGG modules presence/absence...')
            reconstructed_modules <- reconstruct_modules(from = object@annotations_list,
                                                         module_vector = module_vector,
                                                         module_detect_type = module_detect_type)

            reconstruction_results <- return_reconstructed_modules(from = reconstructed_modules,
                                                                   module_detect_type = module_detect_type,
                                                                   output_dir = output_dir,
                                                                   output_prefix = output_prefix,
                                                                   overwrite = overwrite)
            return(reconstruction_results)
          }
)



#' This function is used to reconstruct metabolic pathways based on present KEGG Orthology terms. For any
#' incomplete pathways with missing reactions, it will then calculate the probability that each missing reaction
#' is present in the input genome but was missed in the sampling process. It takes in user input in the form of
#' an object created with a read_data function call.
#' @param object Tibble object created with the read_data function - KEGG orthology annotation data for one or more bacterial genome annotation files.
#' @param module_vector Character vector. An optional character vector of specific KEGG Modules to scan user annotations for, and optionally to predict for. Default is a character vector containing identifiers for all available modules.
#' @param module_detect_type Character. Method 'extract' will give detailed output component containing all extracted KO gene identifiers associated with KEGG modules. Method 'detect' will only detect presence/absence of modules and will not provide detailed output in results.
#' @param output_dir Character. The full or relative path to an output directory where result and summary output will be saved as TSV flatfiles
#' @param output_prefix Character. Optional string to prefix to MetaPredict output files. Default is NULL, resulting in output files with default names.
#' @param overwrite Logical. If TRUE, output files from running the MetaPredict function that have the same name as existing files in the output directory will overwrite those existing files. Default is FALSE.

setMethod(f = 'metapredict',
          signature = 'MetaPredictInput',
          definition = function(object, module_vector = names(all_models), #predict_models = TRUE,
                                module_detect_type = c('extract', 'detect'),
                                output_dir = NULL, output_prefix = NULL, overwrite = FALSE) {

            stopifnot(is(object, 'MetaPredictInput'))

            cli::cli_h1('Starting MetaPredict')
            module_detect_type <- match.arg(module_detect_type)
            object <- check_data(object, module_vector = module_vector)

            cli::cli_alert_info('Reconstructing KEGG modules presence/absence...')
            reconstructed_modules <- reconstruct_modules(from = object@annotations_list,
                                                         module_vector = module_vector,
                                                         module_detect_type = module_detect_type)

            cli::cli_alert_info('Performing prediction calculations...')
            prediction_results <- return_predictions(from = object@annotations_list,
                                                     using = reconstructed_modules,
                                                     module_vector = module_vector,
                                                     module_detect_type = module_detect_type,
                                                     output_dir = output_dir,
                                                     output_prefix = output_prefix,
                                                     overwrite = overwrite)
            return(prediction_results)
          }
)




#setMethod(f = 'save_results.predictions')



#setMethod(f = 'save_results.reconstructions')






# helper functions that return a custom class object ----------------------

return_reconstructed_modules <- function(from,
                                         module_detect_type = module_detect_type,
                                         output_dir = output_dir,
                                         output_prefix = output_prefix,
                                         overwrite = overwrite) {

  summary <- summarize_recon(.recon = from$pres_abs_tbl,
                             .module_metadata = module_metadata)

  if (module_detect_type == 'extract') {
    results <- new('MetaReconstructResult',
                   summary = summary,
                   module_reconstructions = from$pres_abs_tbl,
                   module_extractions = from$ko_extractions)

  } else {
    results <- new('MetaReconstructResult',
                   summary = summary,
                   module_reconstructions = from$pres_abs_tbl)
  }

  if (!is.null(output_dir)) {
    cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
    save_recon(results,
               output_dir,
               output_prefix = output_prefix,
               overwrite = overwrite)
  } else {
    cli::cli_alert_success('All done.')
    cli::cli_alert_info('To save results, use save_recon().')
  }

  return(results)
}



return_predictions <- function(from, using,
                               module_vector = module_vector, module_detect_type = module_detect_type,
                               output_dir = output_dir, output_prefix = output_prefix, overwrite = overwrite) {

  from <- create_kegg_matrix(from) %>%
    predict(caret::preProcess(., method = c('center', 'scale')), .) %>%
    suppressWarnings()

  predictions <- purrr::map2_dfc(all_models[module_vector], names(all_models[module_vector]), ~
                                   named_predict(.y, .x, from)) %>%
    suppressWarnings() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.x))) %>%
    dplyr::mutate(genome_name = using$pres_abs_tbl$genome_name, .before = 1)

  predictions <- put_na(using$pres_abs_tbl, predictions)

  summary <- summarize_results(using$pres_abs_tbl, .pred = predictions, .module_metadata = module_metadata)

  if (module_detect_type == 'extract') {
    results <- new('MetaPredictResult',
                   summary = summary,
                   module_reconstructions = using$pres_abs_tbl,
                   module_predictions = predictions,
                   module_extractions = using$ko_extractions)
  } else {
    results <- new('MetaPredictResult',
                   summary = summary,
                   module_reconstructions = using$pres_abs_tbl,
                   module_predictions = predictions)
  }

  if (!is.null(output_dir)) {
    cli::cli_alert_info('All done. Saving results to directory: {output_dir}')
    save_results(results, output_dir, output_prefix = output_prefix, overwrite = overwrite)
  } else {
    cli::cli_alert_success('All done.')
    cli::cli_alert_info('To save results, use save_results().')
  }

  return(results)
}
