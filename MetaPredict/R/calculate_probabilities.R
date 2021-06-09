#' @importFrom magrittr "%>%"

# p(x_j|k) = CONDITIONAL PROBABILITY that X_j = x_j given that K_j = n
# n = length of module (# genes)
# x_j = THE NUMBER of genes from the module present in an input FRAGMENT f_j (the OBSERVED value of X_j)
# k_j = OBSERVED VALUE OF K_j, which is the NUMBER of MP genes in a genome g_j, of which f_j is a fragment
get_p.xj_n <- function(n, x_j, p_j) {
  p.xj_n <- choose(n, x_j) * p_j^x_j * (1 - p_j)^(n - x_j)
  return(p.xj_n)
}



get_p.xj_k <- function(n, x_j, p_j) {
  k <- x_j:(n - 1)
  p.xj_k <- choose(k, x_j) * p_j^x_j * (1 - p_j)^(k - x_j)
  return(p.xj_k)
}



# p(x_j|nbar) = CONDITIONAL PROBABILITY that X_j = x_j given that K_j != n
#requires
get_p.xj_nbar <- function(n, x_j, p_j, pj.k) {
  p.xj_k <- get_p.xj_k(n, x_j, p_j)
  p.xj_nbar = sum(p.xj_k * pj.k)
  return(p.xj_nbar)
}



## m_l = length of collection C_l
## m_j = length of collection C_j
LogL.betabinomial <- function(params, y_l, n, m_l) {
  alpha <- params[1]
  beta <- params[2]
  sum(-lgamma(alpha + beta) - lgamma(alpha + y_l)
      - lgamma(beta + (n * m_l) - y_l) + lgamma(alpha)
      + lgamma(beta) + lgamma(alpha + beta + (n * m_l)))
}



## m_l = length of collection C_l
## m_j = length of collection C_j
## C_j = collection of known genomes to which fragments f_j belong
## n = length of module (# genes)
## alpha + beta MLEs
## k_j = OBSERVED VALUE OF K_j, which is the NUMBER of MP genes in a genome g_j, of which f_j is a fragment
get_pj.k_betabinomial <- function(alpha, beta, y_j, n, m_j, x_j) { # k = seq_along(x_j:(n - 1)) # x_j should never equal n; if it does, the module step is complete/present
  k <- x_j:(n - 1)
  A_j <- alpha + y_j
  B_j <- beta + (n * m_j) - y_j
  pj.k <- exp(lgamma(A_j + B_j) + lgamma(A_j + k) + lgamma(B_j + n - k) - lgamma(A_j) - lgamma(B_j) - lgamma(A_j + B_j + n))
  return(pj.k)
}



get_pj.n_betabinomial <- function(alpha, beta, y_j, n, m_j, x_j) {
  k <- n
  A_j <- alpha + y_j
  B_j <- beta + (n * m_j) - y_j
  pj.n <- exp(lgamma(A_j + B_j) + lgamma(A_j + k) + lgamma(B_j + n - k) - lgamma(A_j) - lgamma(B_j) - lgamma(A_j + B_j + n))
  return(pj.n)
}



calculate_posterior <- function(p.xj_n, p.xj_nbar, alpha, beta, y_j, n, m_j, x_j) {
  pj.n <- get_pj.n_betabinomial(alpha, beta, y_j, n, m_j, x_j)
  p.n_xj <- (p.xj_n * pj.n) / ((p.xj_n * pj.n) + (p.xj_nbar * (1 - pj.n)))
  return(p.n_xj)
}



# Calculate the posterior probability that a module step is present,
# given the genome completeness estimate p
#' @export
get_posteriors <- function(modules_list, strict = strict) {
  posterior_results <- list()
  if (strict == TRUE) {
    modules_list <- modules_list %>%
      purrr::keep(~ nrow(.x) > 0) %>%
      purrr::set_names(nm = purrr::map_chr(., ~ unique(.x$taxonomy)))
  }
  posterior_results <- purrr::map(seq_along(modules_list), ~ {
    posterior_results[[.x]] <- modules_list[[.x]] %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pj.k = list(get_pj.k_betabinomial(alpha, beta, y_j, n, m_j, x_j)),
             p.xj_nbar = get_p.xj_nbar(n, x_j, p_j, pj.k),
             p.xj_n = list(get_p.xj_n(n, x_j, p_j)),
             probability = calculate_posterior(p.xj_n, p.xj_nbar, alpha, beta, y_j, n, m_j, x_j)) %>%
      dplyr::select(-c(n, x_j, p_j, y_j, m_j, alpha, beta, pj.k, p.xj_nbar, p.xj_n))
  }) %>%
    purrr::set_names(nm = names(modules_list))
  return(posterior_results)
}
