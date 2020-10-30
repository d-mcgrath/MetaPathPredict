#' @importFrom magrittr "%>%"
#' @import dplyr
LogL.betabinomial <- function(v, y_k = 5, n_k = 15) {
  alpha <- v[1]
  beta <- v[2]
  sum(-lgamma(alpha + y_k) - lgamma(beta + n_k - y_k) + lgamma(alpha + beta + n_k)
      + lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta))
}



calculate_p <- function(steps, rxn_matrix, taxonomic_lvl, organism, scan_list, res_list) {
  rlang::enquo(taxonomic_lvl)
  options(dplyr.summarise.inform = FALSE)

  n_k <- (rxn_matrix %>%
            group_by(.data[[!!taxonomic_lvl]]) %>%
            dplyr::filter(.data[[!!taxonomic_lvl]] != organism) %>%
            ungroup() %>%
            group_by(Genus) %>%
            summarize(n_k = length(Genus)))$n_k

  n_j <- (rxn_matrix %>%
            group_by(.data[[!!taxonomic_lvl]]) %>%
            dplyr::filter(.data[[!!taxonomic_lvl]] == organism) %>%
            summarize(n_j = length(.data[[!!taxonomic_lvl]])))$n_j

  y_k <- rxn_matrix %>%
    group_by(.data[[!!taxonomic_lvl]]) %>%
    dplyr::filter(.data[[!!taxonomic_lvl]] != organism) %>%
    ungroup() %>%
    group_by(Genus) %>%
    summarize(across(all_of(steps), sum)) %>%
    select(-Genus)

  y_j <- rxn_matrix %>%
    group_by(.data[[!!taxonomic_lvl]]) %>%
    dplyr::filter(.data[[!!taxonomic_lvl]] == organism) %>%
    summarize(across(all_of(steps), sum)) %>%
    select(-.data[[!!taxonomic_lvl]])

  #if (length(n_k) >= 2) { #MLE using parallelize-able furrr loops
  #  predictions <- furrr::future_map2(1:length(y_k), 1:length(y_j), .progress = TRUE, ~ {
  #    opt <- optim(par = c(0.01, 0.01), LogL.betabinomial, #calculate maximum likelihood estimates of alpha and beta
  #                 y_k = y_k[[.x]],
  #                 n_k = n_k,
  #                 method = 'L-BFGS-B', lower = 1e-10, upper = 1e6)
  #    res <- (opt$par[1] + y_j[[.y]]) / (opt$par[1] + opt$par[2] + n_j)

  #    names(res) <- colnames(y_k)[.x] #y_k and y_j have same colnames
  #    return(res)
  #    })

  #} else {
  #  predictions <- furrr::future_map(1:length(y_j), .progress = TRUE, ~ {
  #    res <- (1 + y_j[[.x]]) / (1 + 1 + n_j) #if >= 2 collections cannot be used to estimate alpha & beta,

  #    names(res) <- colnames(y_j)[.x]
  #    return(res)
  #    })
  #}

  #scan_list <- scan_list %>%
  #  left_join(tibble(reaction = map_chr(predictions, names),
  #                   probability = as_vector(predictions)), by = 'reaction') ###  end of furr loop style







  if (length(n_k) >= 2) { ###MLE estimation using vectorized functions
    estimates <- y_k %>%
      summarize(across(everything(), ~ optim(par = c(0.01, 0.01), LogL.betabinomial, #calculate maximum likelihood estimates of alpha and beta
                                             y_k = .x, n_k = n_k,
                                             method = 'L-BFGS-B', lower = 1e-10, upper = 1e6)))
    estimates <- estimates %>%
      summarize(across(everything(), ~ (.x[[1]][1] + y_j[[cur_column()]]) / (.x[[1]][1] + .x[[1]][2] + n_j) ))

  } else {
    estimates <- y_j %>%
      summarize(across(everything(), ~ (1 + .x) / (1 + 1 + n_j) ))
  }

  scan_list <- scan_list %>%
    select(-probability) %>%
    left_join(estimates %>%
                tidyr::pivot_longer(everything(), names_to = 'step',
                                                values_to = 'probability') %>%
                dplyr::mutate(probability = as.character(probability)),
              by = 'step') #### end of vectorized approach


  res_list <- res_list %>%
    dplyr::filter(!is.na(probability)) %>%
    full_join(scan_list, by = c('organism', 'name', 'module', 'full', 'step', 'probability')) %>%
    select(organism, name, module, full, step, probability) %>%
    arrange(module, step)

  return(res_list)
}



no_optim <- function(steps, rxn_matrix, organism, scan_list, res_list) {
  options(dplyr.summarise.inform = FALSE)

  n_j <- (rxn_matrix %>%
            group_by(Domain) %>%
            summarize(n_j = length(Domain)))$n_j

  y_j <- rxn_matrix %>%
    group_by(Domain) %>%
    summarize(across(all_of(steps), sum)) %>%
    select(-Domain)




  #predictions <- furrr::future_map(1:length(y_j), .progress = T, ~ { ### MLE estimation with parallelizable-processing
  #  res <- (1 + y_j[[.x]]) / (1 + 1 + n_j)

  #  names(res) <- colnames(y_j)[.x]
  #  return(res)
  #  })

  #scan_list <- scan_list %>%
  #  left_join(tibble(reaction = map_chr(predictions, names),
  #                   probability = as_vector(predictions)), by = 'reaction') ###end of furrr loop style





  estimates <- y_j %>% ### MLE with vectorized summarize
    summarize(across(everything(), ~ (1 + .x) / (1 + 1 + n_j) ))


  scan_list <- scan_list %>%
    select(-probability) %>%
    left_join(estimates %>%
                tidyr::pivot_longer(everything(), names_to = 'step',
                                    values_to = 'probability') %>%
                dplyr::mutate(probability = as.character(probability)),
              by = 'step') #### end of vectorized approach


  res_list <- res_list %>%
    dplyr::filter(!is.na(probability)) %>%
    full_join(scan_list, by = c('organism', 'name', 'module', 'full', 'step', 'probability')) %>%
    select(organism, name, module, full, step, probability) %>%
    arrange(module, step)

  return(res_list)
}



taxonomy_not_found <- function(organism, scan_list, res_list) {
  scan_list <- scan_list %>%
    mutate(probability = 'taxonomy not found; no calculation done')

  res_list <- res_list %>%
    dplyr::filter(!is.na(probability)) %>%
    full_join(scan_list, by = c('organism', 'name', 'module', 'full', 'step', 'probability')) %>%
    select(organism, name, module, full, step, probability) %>%
    arrange(module, step)
}



pull_data <- function(steps, organism, scan_list, res_list)  {
  if (organism %in% bacteria.rxn.matrix$Genus) {
    res <- calculate_p(steps, bacteria.rxn.matrix, 'Genus', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Family) {
    res <- calculate_p(steps, bacteria.rxn.matrix, 'Family', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Order) {
    res <- calculate_p(steps, bacteria.rxn.matrix, 'Order', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Class) {
    res <- calculate_p(steps, bacteria.rxn.matrix, 'Class', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Phylum) {
    res <- calculate_p(steps, bacteria.rxn.matrix, 'Phylum', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Domain) {
    res <- no_optim(steps, bacteria.rxn.matrix, organism, scan_list, res_list)


  #} else if (organism %in% archaea.rxn.matrix$Genus) {
  #  res <- calculate_p(steps, archaea.rxn.matrix, 'Genus', organism, scan_list, res_list)

  #} else if (organism %in% archaea.rxn.matrix$Family) {
  #  res <- calculate_p(steps, archaea.rxn.matrix, 'Family', organism, scan_list, res_list)

  #} else if (organism %in% archaea.rxn.matrix$Order) {
  #  res <- calculate_p(steps, archaea.rxn.matrix, 'Order', organism, scan_list, res_list)

  #} else if (organism %in% archaea.rxn.matrix$Class) {
  #  res <- calculate_p(steps, archaea.rxn.matrix, 'Class', organism, scan_list, res_list)

  #} else if (organism %in% archaea.rxn.matrix$Phylum) {
  #  res <- calculate_p(steps, archaea.rxn.matrix, 'Phylum', organism, scan_list, res_list)

  #} else if (organism %in% archaea.rxn.matrix$Domain) {
  #  res <- no_optim(steps, archaea.rxn.matrix, scan_list, res_list)

  } else {
    res <- taxonomy_not_found(organism, scan_list, res_list)
  }
  return(res)
}



summarize_output <- function(results) {
  results <- results %>%
    dplyr::filter(!(duplicated(step))) %>%
    tally(probability == 'Present') %>%
    full_join(suppressWarnings(
      results %>%
        group_by(module, step) %>%
        dplyr::filter(!any(probability == 'Present')) %>%
        mutate(probability = as.numeric(probability)) %>%
        dplyr::filter(!is.na(probability)) %>%
        ungroup() %>%
        group_by(module) %>%
        dplyr::filter(!(duplicated(step))) %>%
        tally(length(module))), by = 'module') %>%
    full_join(results %>%
                group_by(module, step) %>%
                dplyr::filter(!any(probability == 'Present')) %>%
                suppressWarnings(mutate(probability = as.numeric(probability))) %>%
                dplyr::filter(!is.na(probability)) %>%
                ungroup() %>%
                group_by(module) %>%
                dplyr::filter(!(duplicated(step))) %>%
                tally(probability >= 0.75), by = 'module') %>%
    full_join(results %>%
                dplyr::filter(!(duplicated(step))) %>%
                tally(length(module)), by = 'module') %>%
    rename(module_steps_present = 2, module_steps_predicted = 3,
           `probability > 75%` = 4, total_module_steps = 5) %>%
    mutate(module_steps_predicted = case_when(
      is.na(module_steps_predicted) ~ 0L,
      TRUE ~ module_steps_predicted),
      `probability > 75%` = case_when(
        is.na(`probability > 75%`) ~ 0L,
        TRUE ~ `probability > 75%`)) %>%
    inner_join(distinct(select(results, module, name)), by = 'module') %>%
    select(module, name, module_steps_present,
           module_steps_predicted, `probability > 75%`, total_module_steps) %>%
    mutate(percent_present =
             paste0(signif(module_steps_present/total_module_steps, digits = 3) * 100, '%'),
           `predicted_completeness (probability>75%)` =
             paste0(signif((`probability > 75%` +
                              module_steps_present) /
                             total_module_steps, digits = 3) * 100, '%'),
           module_steps_present =
             paste0(module_steps_present, '/', total_module_steps),
           module_steps_predicted =
             paste0(module_steps_predicted, '/', total_module_steps),
           `probability > 75%` =
             paste0(`probability > 75%`, '/', total_module_steps))
  return(results)
}
