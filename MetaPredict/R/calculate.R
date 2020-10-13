#' @importFrom magrittr "%>%"
#' @import dplyr
LogL.betabinomial <- function(v, y_k = 5, n_k = 15) {
  alpha <- v[1]
  beta <- v[2]
  sum(-lgamma(alpha + y_k) - lgamma(beta + n_k - y_k) + lgamma(alpha + beta + n_k)
    + lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta))
}



calculate_p <- function(reactions, rxn_matrix, taxonomic_lvl, organism, scan_list, res_list) {
  rlang::enquo(taxonomic_lvl)
  options(dplyr.summarise.inform = FALSE)

  collection.k <- rxn_matrix %>%
    group_by(.data[[!!taxonomic_lvl]]) %>%
    filter(.data[[!!taxonomic_lvl]] != organism) %>%
    ungroup() %>%
    group_by(Genus)

  collection.j <- rxn_matrix %>%
    group_by(.data[[!!taxonomic_lvl]]) %>%
    filter(.data[[!!taxonomic_lvl]] == organism)

  n_k <- (collection.k %>%
            summarize(n_k = length(Genus)))$n_k

  n_j <- (collection.j %>%
            summarize(n_j = length(.data[[!!taxonomic_lvl]])))$n_j

  y_k <- collection.k %>%
    summarize(across(all_of(reactions), sum)) %>%
    select(-Genus)

  y_j <- collection.j %>%
    summarize(across(all_of(reactions), sum)) %>%
    select(-.data[[!!taxonomic_lvl]])


  predictions <- furrr::future_map2(1:length(y_k), 1:length(y_j), .progress = T, ~ {
    if (length(n_k) >= 2) {
      opt <- optim(par = c(0.01, 0.01), LogL.betabinomial,
                   y_k = y_k[[.x]],
                   n_k = n_k,
                   method = 'L-BFGS-B', lower = 1e-10, upper = 1e6)
      res <- (opt$par[1] + y_j[[.y]]) / (opt$par[1] + opt$par[2] + n_j)

    } else {
      res <- NA }

    names(res) <- colnames(y_k)[.x] #y_k and y_j have same colnames
    return(res)
  })

  scan_list <- scan_list %>%
    left_join(tibble(reaction = map_chr(predictions, names),
                     probability = as_vector(predictions)), by = 'reaction')

  res_list <- res_list %>%
    full_join(scan_list, by = c('ko_description', 'reaction', 'reaction_description',
                                'pathway', 'pathway_name',
                                'pathway_class', 'organism')) %>%
    mutate(probability = as.character(probability), probability = case_when(
      !(is.na(ko_term)) & is.na(probability) ~ 'Present',
      TRUE ~ probability), organism = case_when(
        is.na(organism) ~ unique(na.omit(organism)),
        TRUE ~ organism)) %>%
    select(organism, pathway, reaction, ko_term, probability,
           pathway_name, reaction_description, ko_description, pathway_class)

  message('Done with ', unique(na.omit(res_list$organism)))
  return(res_list)
}



no_optim <- function(reactions, rxn_matrix) {
  options(dplyr.summarise.inform = FALSE)

  collection.j <- rxn_matrix %>%
    group_by(Genus)

  n_j <- (collection.j %>%
            summarize(n_j = length(Genus)))$n_j

  y_j <- collection.j %>%
    summarize(across(all_of(reactions), sum)) %>%
    select(-Genus)

  predictions <- furrr::future_map(1:length(y_j), .progress = T, ~ {
    res <- (1 + y_j[[.x]]) / (1 + 1 + n_j)

    names(res) <- colnames(y_j)[.x]
    return(res)
    })

  scan_list <- scan_list %>%
    left_join(tibble(reaction = map_chr(predictions, names),
                     probability = as_vector(predictions)), by = 'reaction')

  res_list <- res_list %>%
    full_join(scan_list, by = c('ko_description', 'reaction', 'reaction_description',
                                'pathway', 'pathway_name',
                                'pathway_class', 'organism')) %>%
    mutate(probability = as.character(probability), probability = case_when(
      !(is.na(ko_term)) & is.na(probability) ~ 'Present',
      TRUE ~ probability), organism = case_when(
        is.na(organism) ~ unique(na.omit(organism)),
        TRUE ~ organism)) %>%
    select(organism, pathway, reaction, ko_term, probability,
           pathway_name, reaction_description, ko_description, pathway_class)

  message('Done with ', unique(na.omit(res_list$organism)), '\n')
  return(res_list)
}



pull_data <- function(reactions, organism, scan_list, res_list)  {
  if (organism %in% bacteria.rxn.matrix$Genus) {
    res <- calculate_p(reactions, bacteria.rxn.matrix, 'Genus', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Family) {
    res <- calculate_p(reactions, bacteria.rxn.matrix, 'Family', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Order) {
    res <- calculate_p(reactions, bacteria.rxn.matrix, 'Order', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Class) {
    res <- calculate_p(reactions, bacteria.rxn.matrix, 'Class', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Phylum) {
    res <- calculate_p(reactions, bacteria.rxn.matrix, 'Phylum', organism, scan_list, res_list)

  } else if (organism %in% bacteria.rxn.matrix$Domain) {
    res <- no_optim(reactions, bacteria.rxn.matrix, scan_list, res_list)


  } else if (organism %in% imgm.archaea.rxn.matrix$Genus) {
    res <- calculate_p(reactions, imgm.archaea.rxn.matrix, 'Genus', organism, scan_list, res_list)

  } else if (organism %in% imgm.archaea.rxn.matrix$Family) {
    res <- calculate_p(reactions, imgm.archaea.rxn.matrix, 'Family', organism, scan_list, res_list)

  } else if (organism %in% imgm.archaea.rxn.matrix$Order) {
    res <- calculate_p(reactions, imgm.archaea.rxn.matrix, 'Order', organism, scan_list, res_list)

  } else if (organism %in% imgm.archaea.rxn.matrix$Class) {
    res <- calculate_p(reactions, imgm.archaea.rxn.matrix, 'Class', organism, scan_list, res_list)

  } else if (organism %in% imgm.archaea.rxn.matrix$Phylum) {
    res <- calculate_p(reactions, imgm.archaea.rxn.matrix, 'Phylum', organism, scan_list, res_list)

  } else if (organism %in% imgm.archaea.rxn.matrix$Domain) {
    res <- no_optim(reactions, imgm.archaea.rxn.matrix, scan_list, res_list)

  } else {
    stop('Organism taxonomy not detected for ', organism, ', even at the domain level. Please see MetaPredict usage instructions [-h].') }

  return(res)
}
