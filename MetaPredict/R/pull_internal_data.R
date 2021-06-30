#' @importFrom magrittr "%>%"
#' @export
get_gene_counts_from <- function(userData) {
  input_data <- userData %>%
    dplyr::ungroup() %>%
    {if ('genome_name' %in% colnames(.))
      dplyr::select(., k_number, gene_name, genome_name)
      else if ('metagenome_name' %in% colnames(.))
        dplyr::select(., k_number, gene_name, taxonomy)
      else stop(cli::cli_alert_danger('Error: either "genome_name" or "metagenome_name" column is missing from the input data.'))
    } %>%
    dplyr::rowwise() %>%
    dplyr::mutate(k_numbers = list(stringr::str_extract_all(
      string = k_number,
      pattern = patt.kegg_modules$k_numbers)),
      k_numbers = list(purrr::map(k_numbers, ~
                                    if (length(.x) == 0) {.x <- NA_character_
                                    } else {paste0(.x, collapse = ';')})),
      x_j.k_numbers = list(dplyr::tibble(k_numbers) %>%
                             dplyr::mutate(step = patt.kegg_modules$step) %>%
                             dplyr::group_by(step) %>%
                             dplyr::summarize(k_numbers = paste0(k_numbers[!is.na(k_numbers)], collapse = ' '), .groups = 'drop') %>%
                             dplyr::select(-step) %>%
                             dplyr::mutate(k_numbers = stringr::str_replace(k_numbers, stringr::regex('^$'), NA_character_),
                                           x_j = dplyr::if_else(is.na(k_numbers), 0, (stringr::str_count(k_numbers, pattern = ' ') + 1))))) %>%
    dplyr::ungroup()

  gene_counts.list <- list()
  gene_counts.list <- purrr::map(1:length(input_data$x_j.k_numbers), ~ {
    gene_counts.list[[.x]] <- input_data$x_j.k_numbers[[.x]]
  }) %>%
    {if ('genome_name' %in% colnames(userData))
      purrr::set_names(., nm = userData$genome_name)
      else if ('metagenome_name' %in% colnames(userData))
        purrr::set_names(., nm = input_data$taxonomy)
      else stop(cli::cli_alert_danger("Error: 'genome_name' or 'metagenome_name' column is missing from input data."))
    }
  return(gene_counts.list)
}


dt.get_gene_counts_from <- function (userData) {
  userData <- lazy_dt(userData) # this is NOT done yet. want to rewrite this w/ dtplyr or data.table

  input_data <- userData %>%
    dplyr::ungroup() %>%
    {if ('genome_name' %in% colnames(.))
      dplyr::select(., k_number, gene_name, genome_name)
      else if ('metagenome_name' %in% colnames(.))
        dplyr::select(., k_number, gene_name, taxonomy)
      else stop(cli::cli_alert_danger('Error: either "genome_name" or "metagenome_name" column is missing from the input data.'))
    } %>%
    dplyr::rowwise() %>%
    dplyr::mutate(k_numbers = list(stringr::str_extract_all(
      string = k_number,
      pattern = patt.kegg_modules$k_numbers)),
      k_numbers = list(purrr::map(k_numbers, ~
                                    if (length(.x) == 0) {.x <- NA_character_
                                    } else {paste0(.x, collapse = ';')})),
      x_j.k_numbers = list(dplyr::tibble(k_numbers) %>%
                             dplyr::mutate(step = patt.kegg_modules$step) %>%
                             dplyr::group_by(step) %>%
                             dplyr::summarize(k_numbers = paste0(k_numbers[!is.na(k_numbers)], collapse = ' '), .groups = 'drop') %>%
                             dplyr::select(-step) %>%
                             dplyr::mutate(k_numbers = stringr::str_replace(k_numbers, stringr::regex('^$'), NA_character_),
                                           x_j = dplyr::if_else(is.na(k_numbers), 0, (stringr::str_count(k_numbers, pattern = ' ') + 1))))) %>%
    dplyr::ungroup()

  gene_counts.list <- list()
  gene_counts.list <- purrr::map(1:length(input_data$x_j.k_numbers), ~ {
    gene_counts.list[[.x]] <- input_data$x_j.k_numbers[[.x]]
  }) %>%
    {if ('genome_name' %in% colnames(userData))
      purrr::set_names(., nm = userData$genome_name)
      else if ('metagenome_name' %in% colnames(userData))
        purrr::set_names(., nm = input_data$taxonomy)
      else stop(cli::cli_alert_danger("Error: 'genome_name' or 'metagenome_name' column is missing from input data."))
    }
  return(gene_counts.list)
}


# previous version, deprecated: 04/17/21
#get.x_j <- function(userData) {
#  input_data <- userData %>%
#    dplyr::ungroup() %>%
#    {if ('genome_name' %in% colnames(userData))
#      dplyr::select(., k_number, genome_name)
#      else dplyr::select(., k_number, taxonomy)
#    } %>%
#    dplyr::rowwise() %>%
#    dplyr::mutate(x_j = list(as.integer(stringr::str_detect(
#      string = k_number,
#      pattern = patt.kegg_modules$k_numbers)))) %>%
#    dplyr::ungroup() %>%
#    dplyr::mutate(x_j = purrr::map(x_j, ~
#                                     dplyr::tibble(x_j = .x) %>%
#                                     dplyr::mutate(step = patt.kegg_modules$step) %>%
#                                     dplyr::group_by(step) %>%
#                                     dplyr::summarize(x_j = sum(x_j), .groups = 'drop') %>%
#                                     dplyr::select(-step)))
#  gene_counts.list <- list()
#  gene_counts.list <- purrr::map(1:length(input_data$x_j), ~ {
#    gene_counts.list[[.x]] <- input_data$x_j[[.x]][[1]]
#  }) %>%
#    {if ('genome_name' %in% colnames(userData))
#      purrr::set_names(., nm = userData$genome_name)
#      else if ('metagenome_name' %in% colnames(userData))
#    purrr::set_names(., nm = input_data$taxonomy)
#      else stop(cli::cli_alert_danger("Error: 'genome_name' or 'metagenome_name' column is missing from input data."))
#    }
#  return(gene_counts.list)
#}



#previous version - deprecated
#get.x_j <- function(userData) {
#  input_data <- userData %>%
#    dplyr::ungroup() %>%
#    {if ('genome_name' %in% colnames(userData))
#      dplyr::select(., k_number, genome_name) %>%
#        tibble::column_to_rownames(var = 'genome_name')
#      else dplyr::select(., k_number, taxonomy) %>%
#        tibble::column_to_rownames(var = 'taxonomy')} %>%
#    t() %>%
#    dplyr::as_tibble(~ vctrs::vec_as_names(repair = "unique", quiet = TRUE))
#
#  x_j.tibble <- patt.kegg_modules %>%
#    dplyr::bind_cols(input_data) %>% # binding COLUMNS, not ROWS of the input data with the patt.kegg_modules tibble
#    dplyr::group_by(step, k_numbers) %>%
#    dplyr::summarize(dplyr::across(!1:3, ~ as.integer(stringr::str_detect(
#      string = .x, pattern = k_numbers))), .groups = 'keep') %>%
#    dplyr::ungroup(k_numbers) %>%
#    dplyr::select(-c(k_numbers)) %>%
#    dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x)), .groups = 'drop') %>%
#    dplyr::select(-step)
#
#  gene_counts.list <- list()
#  gene_counts.list <- purrr::map(1:length(x_j.tibble), ~ {
#    gene_counts.list[[.x]] <- x_j.tibble[[.x]]
#  }) %>%
#    purrr::set_names(nm = names(x_j.tibble))
#  return(gene_counts.list)
#}



#' @export
map_modules_to <- function(gene_counts.list, userData, strict = FALSE) {
  modules_list <- list()
  modules_list <- purrr::map(1:length(gene_counts.list), ~ {
    modules_list[[.x]] <- all_kegg_modules %>%
      dplyr::ungroup() %>%
      #dplyr::mutate(x_j = gene_counts.list[[.x]]) %>%
      dplyr::bind_cols(gene_counts.list[[.x]]) %>%
      dplyr::relocate(k_numbers, .after = rule) %>%
      dplyr::mutate(taxonomy = userData$taxonomy[.x], .before = 1) %>%
      {if (unique(userData$data_type) == 'genome') dplyr::mutate(., 'genome_name' = userData$genome_name[.x], .after = 1)
        else if (unique(userData$data_type) == 'metagenome') dplyr::mutate(., 'metagenome_name' = userData$metagenome_name[.x], .after = 1)
        else stop (cli::cli_alert_danger('Error: data_type was not interpreted; input data must be one or more genomes or metagenomes.'))} %>%
      dplyr::mutate(domain = get_domain(userData$taxonomy[.x]), .after = 1) %>%
      get_lowest_taxonomy() %>%
      dplyr::mutate(module_step_present = dplyr::case_when(n == x_j ~ TRUE, n != x_j ~ FALSE),
                    p_j = userData$completeness[.x]) %>%
      dplyr::arrange(step) %>%
      dplyr::group_by(module) %>%
      {if (strict == TRUE) dplyr::filter(., !all(module_step_present == FALSE)) else (.)} # filters out pathways for which there is no evidence of their presence in the gene annotations
  }) %>%
    purrr::set_names(nm = names(gene_counts.list))
  return(modules_list)
}



get_domain <- function(taxonomy) { #NOTE: some columns of bact_ & arch_phylo_key contain NA, hence the na.omit()
  if (purrr::some(1:7, ~ {
    any(stringr::str_detect(na.omit(bact_phylo_key[[.x]]), stringr::regex(paste0('^', taxonomy, '$'), ignore_case = TRUE)))})) {
    domain = 'Bacteria'
  } else if (purrr::some(1:7, ~ { #change to 1:7 when arch_phylo_key gets species-level support
    any(stringr::str_detect(na.omit(arch_phylo_key[[.x]]), stringr::regex(paste0('^', taxonomy, '$'), ignore_case = TRUE)))})) {
    domain = 'Archaea'
  } else {
    domain = 'Unknown'
  }
  return(domain)
}



get_lowest_taxonomy <- function(.data) {
  if (unique(.data$domain) == 'Bacteria') {
    search_tree(.data, bact_phylo_key, .list_num = 1, .domain = 'Bacteria_')
  } else if (unique(.data$domain) == 'Archaea') {
    search_tree(.data, arch_phylo_key, .list_num = 2, .domain = 'Archaea_')
  } else {
    dplyr::mutate(.data, taxonomy_used = 'Unknown')
  }
}



search_tree <- function(.data, taxonomy_key, .list_num, .domain) {
  indices <- taxonomy_key %>% #find out exactly where taxonomy matches are in taxonomic table keys
    purrr::map(~ {stringr::str_detect(.x, stringr::regex(paste0('^', unique(.data$taxonomy), '$'), ignore_case = TRUE))}) %>%
    purrr::map(which) %>%
    purrr::keep(~ !is_empty(.x))

  if (!purrr::is_empty(indices)) { #extract taxonomy from key, including higher levels up to root (domain)
    possible_taxa <- taxonomy_key %>%
      dplyr::slice(purrr::pluck(indices, 1, 1)) %>%
      dplyr::select(purrr::pluck(indices, names, 1):7)

    for (.x in seq_along(possible_taxa)) {
      if (is.na(possible_taxa[[.x]][1])) {
        next
      } else {
        if (any(stringr::str_detect(
          colnames(priors_list[[.list_num]][[paste0(.domain, names(possible_taxa)[.x])]]),
          possible_taxa[[.x]][1]))) {
          .data <- .data %>%
            dplyr::mutate(taxonomy_used = possible_taxa[[.x]][1])
          break
        } else {
          .data <- .data %>%
            dplyr::mutate(taxonomy_used = 'Unknown')
        }
      }
    }
  } else {
    .data <- .data %>%
      dplyr::mutate(taxonomy_used = 'Unknown')
  }
  return(.data)
}



## pulls y_j, m_j, alpha, and beta
#' @export
get_parameters <- function(data_list, moduleVector = NULL, strict = strict)  {
  if (strict == TRUE) {
    data_list <- purrr::keep(data_list, ~ nrow(.x) > 0)
  }
  input_list <- list(bact_list = data_list %>% purrr::keep(~ .x$domain[1] == 'Bacteria'),
                     arch_list = data_list %>% purrr::keep(~ .x$domain[1] == 'Archaea'),
                     unk_list = data_list %>% purrr::keep(~ .x$domain[1] == 'Unknown'))
  res_list <- list() #1:2 for bact and arch_lists; unk needs its own section

  purrr::map_if(seq_along(input_list), ~ length(input_list[[.x]]) != 0, ~ {
    res_list[[.x]] <- list()

    purrr::map(1:length(input_list[[.x]]), function(cur) {
      cur_org <- input_list[[.x]][[cur]]
      cur_org <- cur_org %>% #cur_org is a grouped tibble for the current MAG/metagenome contig cluster
        dplyr::ungroup() #%>%
      #{if (!is.null(moduleVector)) dplyr::filter(., module %in% moduleVector) else (.)}

      cur_org.list <- list(present = dplyr::filter(cur_org, module_step_present == TRUE),
                           absent = dplyr::filter(cur_org, module_step_present == FALSE)) # split cur_org into a list of 2 tibbles based on module presence/absence

      #message('Working on: ', cur_org$lowest[1])

      full_indices <- purrr::map_depth(priors_list[[.x]], .depth = 1, function(.cur_priors) {
        purrr::map_lgl(colnames(.cur_priors), function(.column_name) {
          stringr::str_detect(.column_name, stringr::regex(paste0('^', cur_org$taxonomy_used[1], '$')))
        }) %>%
          which()
      })
      #1st: gets position of list element with first regex match; 2nd: gets position of column within list element containing first regex match
      lowest_index <- c(purrr::detect_index(full_indices, ~ !purrr::is_empty(.x)), purrr::detect(full_indices, ~ !purrr::is_empty(.x)))

      cur_org.list[['absent']] <- cur_org.list[['absent']] %>%
        dplyr::mutate(y_j = get_yj(modules_contained_y[[names(full_indices)[lowest_index[1]]]], cur_org$taxonomy_used[1], step_col = .$step),
                      m_j = get_mj(collection_lengths_m[[names(full_indices)[lowest_index[1]]]], cur_org$taxonomy_used[1])) %>%
        dplyr::left_join(dplyr::select(priors_list[[.x]][[names(full_indices)[lowest_index[1]]]], step, lowest_index[2]) %>%
                           tidyr::separate_rows(dplyr::everything(), sep = ' '), by = 'step') %>%
        #dplyr::rename(parameters = 16) %>% ### POSITION-BASED STEP- MUST CHANGE IF CODE IS UPDATED
        dplyr::rename(parameters = 17) %>%
        tidyr::separate(parameters, into = c('alpha', 'beta'), sep = ';') %>%
        dplyr::mutate(dplyr::across(c(alpha, beta), ~ as.numeric(.x)))

      cur_org.list[['present']] <- cur_org.list[['present']] %>%
        dplyr::mutate(y_j = NA, m_j = NA, alpha = NA, beta = NA)

      res_list[[.x]][[cur]] <- cur_org.list[['absent']] %>%
        dplyr::bind_rows(cur_org.list[['present']]) %>%
        dplyr::arrange(module, step) %>%
        {if (!is.null(moduleVector)) dplyr::filter(., stringr::str_detect(module, stringr::regex(paste0('^', moduleVector, '$', collapse = '|'))))
          else (.)}
        #{if (!is.null(moduleVector)) dplyr::filter(., module %in% moduleVector) else (.)}

      #return(res_list)
    })
  }) %>%
    rlang::squash() %>%
    suppressWarnings() %>%
    purrr::keep(~ is_tibble(.x)) %>%
    purrr::set_names(nm = names(data_list))
}


get_yj <- function(y_list, taxonomy, step_col = NULL) {
  y_j <- y_list %>%
    dplyr::select(dplyr::matches(paste0('^', taxonomy, '$'), ignore.case = TRUE), step, module) %>%
    dplyr::filter(step %in% step_col) %>%
    dplyr::select(-c(module, step))
  return(y_j[[1]])
}



get_mj <- function(m_list, taxonomy) {
  m_j <- m_list %>%
    dplyr::filter(.[[1]] == taxonomy) %>%
    select(-1)
  return(m_j[[1]])
}
