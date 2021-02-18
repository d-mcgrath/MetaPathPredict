#' @importFrom magrittr "%>%"
get.knumber_vectors <- function(userData, posterior_data) { #posterior_data is a list - with posterior probabilities calculated
  input_data <- userData %>%
    dplyr::ungroup() %>%
    {if ('genome_name' %in% colnames(userData))
      dplyr::select(., k_number, genome_name) %>%
        tibble::column_to_rownames(var = 'genome_name')
      else dplyr::select(., k_number, taxonomy) %>%
        tibble::column_to_rownames(var = 'taxonomy')} %>%
    t() %>%
    dplyr::as_tibble(~ vctrs::vec_as_names(repair = "unique", quiet = TRUE))

  knumber_data <- patt.kegg_modules %>%
    dplyr::bind_cols(input_data) %>%
    dplyr::group_by(step, k_numbers) %>%
    dplyr::summarize(dplyr::across(!1:3, ~ stringr::str_extract(
      string = .x, pattern = k_numbers)), .groups = 'drop') %>%
    dplyr::select(-c(k_numbers, step))

  x_j.tibble_list <- list()
  x_j.tibble_list <- purrr::map(1:length(knumber_data), ~ {
    x_j.tibble_list[[.x]] <- knumber_data[.x] %>%
      dplyr::rename(k_number = 1) %>%
      dplyr::mutate(step = patt.kegg_modules$step,
             k_rule = patt.kegg_modules$k_numbers) %>%
      dplyr::left_join(dplyr::select(posterior_data[[.x]], step, `Module step present`, probability), by = 'step') %>%
      dplyr::mutate(probability = dplyr::case_when(`Module step present` == TRUE & is.na(probability) ~ -1,
                                     TRUE ~ probability),
             k_number = dplyr::case_when(is.na(k_number) ~ k_rule,
                                  TRUE ~ k_number)) %>%
      tidyr::separate_rows(k_number, sep = '\\|')

    x_j.tibble_list[[.x]] <- purrr::as_vector(x_j.tibble_list[[.x]]$probability) %>%
      purrr::set_names(x_j.tibble_list[[.x]]$k_number)
  }) %>%
    purrr::set_names(nm = names(knumber_data))
  return(x_j.tibble_list)
}



#' @export
create_kegg_maps <- function(k_vector.list, pathway_id, output_dir = '.',
                             k_numbers = TRUE, ec_numbers = FALSE,
                             present_gene_color = 'plum2', probability_gene_color = 'seagreen2') {
  pathview_output <- list()

  if (k_numbers == TRUE & ec_numbers == FALSE) {
    purrr::imap(k_vector.list, ~ {
      purrr::map(pathway_id, function(.p_id) {
        pathview_output[[.y]] <- pathview::pathview(gene.data = .x, pathway.id = .p_id,
                                                    species = 'ko', out.suffix = names(k_vector.list[.y]),
                                                    kegg.dir = output_dir, low = list(gene = present_gene_color, cpd = 'blue'),
                                                    high = list(gene = probability_gene_color, cpd = 'yellow'),
                                                    bins = list(genes = 14, cpd = 14), map.symbol = TRUE, same.layer = FALSE)
      })
    })
  } else if (k_numbers == FALSE & ec_numbers == TRUE) {
    purrr::imap(k_vector.list, ~ {
      purrr::map(pathway_id, function(.p_id) {
        pathview_output[[.y]] <- pathview::pathview(gene.data = .x, pathway.id = .p_id,
                                                    species = 'ko', out.suffix = names(k_vector.list[.y]),
                                                    kegg.dir = output_dir, low = list(gene = present_gene_color, cpd = 'blue'),
                                                    high = list(gene = probability_gene_color, cpd = 'yellow'),
                                                    bins = list(genes = 14, cpd = 14))
      })
    })
  } else {cli::cli_alert_danger("Argument 'k_numbers' or 'ec_numbers' must be TRUE. Both cannot be TRUE, and both cannot be FALSE. Please make sure only one of them is TRUE.")
  stop()}
}
