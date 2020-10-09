#load required packages and their dependencies
if(all(c('tidyverse', 'furrr', 'optparse') %in% rownames(installed.packages()))) {
  invisible(lapply(c('tidyverse', 'furrr', 'optparse'), function(package) {
    suppressWarnings(suppressPackageStartupMessages(library(package, character.only = T)))})) 
  } else {
  stop('One or more of the required packages are not installed. Are you currently running the 
       MetaPredict conda environment?')
}


#function to parse in user data for pathway reconstruction & reaction prediction
read_data = function(filePath, filePattern) {
  setwd(filePath)
  index = list.files(path = filePath, pattern = filePattern)
  res = list()
  
  res = future_map(1:length(index), .progress = T, ~ {
    col = sub('(.*)-ko.tsv', '\\1', index[.x], perl = T)
    res[[.x]] = read_delim(index[.x], col_types = cols(), delim = '\t') %>%
      slice(-1) %>%
      select(KO, `E-value`) %>%
      mutate(`E-value` = as.numeric(`E-value`)) %>%
      filter(`E-value` <= argv$`e-value`) %>% #e-value is set by argv, default = 1e-3
      select(KO) %>%
      rename(!!col := KO) %>%
      filter(!(duplicated(.data[[col]]))) })
  
  res = res %>%
    tibble() %>%
    unnest(cols = everything()) %>% 
    pivot_longer(cols = everything(), names_to = 'organism', 
                 values_to = 'ko_term', values_drop_na = T) %>%
    group_by(organism)
  
  return(res)
}



#beta binomial pmf/maximum likelihood function
LogL.bb.5 = function(v, y_k = 5, n_k = 15) { alpha = v[1]; beta = v[2];
sum(-lgamma(alpha + y_k) - lgamma(beta + n_k - y_k) + lgamma(alpha + beta + n_k)
    + lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta)) 
}



calculate_sql = function(reactions, organism)  {
  
  bac.table = tbl(meta.DB, 'bacteria.rxn.table')
  arc.table = tbl(meta.DB, 'archaea.rxn.table')
  
  future_map(reactions, .progress = T, ~ {
    if (organism %in% bac.table$Genus) {
      coll.k = bac.table %>%
        filter(Genus != organism) %>%
        summarize(y_k = sum(.data[[.x]]), n_k = length(Genus)) 
      
      coll.j = bac.table %>%
        filter(Genus == organism) %>%
        summarize(y_j = sum(.data[[.x]]), n_j = length(Genus))
      
    } else if (organism %in% arc.table$Genus) {
      coll.k = arc.table %>%
        filter(Genus != organism) %>%
        summarize(y_k = sum(.data[[.x]]), n_k = length(Genus)) 
      
      coll.j = arc.table %>%
        filter(Genus == organism) %>%
        summarize(y_j = sum(.data[[.x]]), n_j = length(Genus))
      
    } else {
      stop('Genus not found. Please see MetaPredict usage instructions [-h].') }
    
    opt = optim(par = c(0.01, 0.01), LogL.bb.5, y_k = coll.k$y_k, n_k = coll.k$n_k,
                method = 'L-BFGS-B', lower = 1e-10, upper = 1e6)
    
    res = (opt$par[1] + coll.j$y_j) / (opt$par[1] + opt$par[2] + coll.j$n_j)
    names(res) = .x
    return(res) })
}



#updated function to find maximum likelihood estimates of alpha and beta
calculate_p = function(reactions, rxn_matrix, taxonomic_lvl, organism, scan_list, res_list) {
  predictions = future_map(reactions, .progress = T, ~ {
    rlang::enquo(taxonomic_lvl)
    
    collection.k = rxn_matrix %>%
      group_by(.data[[!!taxonomic_lvl]]) %>%
      filter(.data[[!!taxonomic_lvl]] != organism) %>%
      #{{taxonomic_lvl}} != regex('unclassified', ignore_case = T)) %>%
      ungroup() %>%
      group_by(Genus) %>%
      summarize(y_k = sum(.data[[.x]]), n_k = length(Genus))
    
    collection.j = rxn_matrix %>%
      group_by(.data[[!!taxonomic_lvl]]) %>%
      filter(.data[[!!taxonomic_lvl]] == organism) %>%
      summarize(y_j = sum(.data[[.x]]), n_j = length(.data[[!!taxonomic_lvl]]))
    
    if (length(collection.k$n_k) >= 2) {
      opt = optim(par = c(0.01, 0.01), LogL.bb.5, y_k = collection.k$y_k, 
                  n_k = collection.k$n_k, method = 'L-BFGS-B', lower = 1e-10, upper = 1e6)
      
      res = (opt$par[1] + collection.j$y_j) / (opt$par[1] + opt$par[2] + collection.j$n_j)
      
    } else {
      res = NA } #need n_k to be >= 2, result is NA otherwise - can't be calculated
    #stop('Cannot calculate MLE for alpha and beta when number of collections k < 2.') }
    
    names(res) = .x #output is a named vector(verify?) with p_j & rxn name for each element
    
    return(res) })
  
  scan_list = scan_list %>%
    mutate(probability = as_vector(predictions)) %>%
    group_by(pathway)
  
  res_list = res_list %>%
    full_join(scan_list, by = c('reaction', 'reaction_description', 
                                'pathway', 'pathway_name', 
                                'pathway_class', 'organism')) %>%
    mutate(probability = as.character(probability), probability = case_when(
      !(is.na(ko_term)) & is.na(probability) ~ 'Present',
      TRUE ~ probability), organism = case_when(
        is.na(organism) ~ unique(na.omit(organism)),
        TRUE ~ organism))
  
  message('Saving output...')
  write_tsv(res_list, path = paste(argv$output, unique(na.omit(res_list$organism)),
                                   '-MetaPredict.tsv', sep = ''))
  message('Done with ', unique(na.omit(res_list$organism)), '\n')
}



no_optim = function(reactions, rxn_matrix) {
  future_map(reactions, .progress = T, ~ {
    
    collection.j = rxn_matrix %>%
      group_by(Genus) %>%
      summarize(y_j = sum(.data[[.x]]), n_k = length(Genus))
    
    #future_map(reactions, .progress = T, ~ {
    res = (1 + collection.j$y_j) / (1 + 1 + collection.j$n_j)
    names(res) = .x
    return(res) })
}



pull_data = function(reactions, organism, scan_list, res_list)  { #added scan argument; to all below as well
  if (organism %in% bacteria.rxn.matrix$Genus) {
    res = calculate_p(reactions, bacteria.rxn.matrix, 'Genus', organism, scan_list, res_list)
    
  } else if (organism %in% bacteria.rxn.matrix$Family) {
    res = calculate_p(reactions, bacteria.rxn.matrix, 'Family', organism, scan_list, res_list)
    
  } else if (organism %in% bacteria.rxn.matrix$Order) {
    res = calculate_p(reactions, bacteria.rxn.matrix, 'Order', organism, scan_list, res_list)
    
  } else if (organism %in% bacteria.rxn.matrix$Class) {
    res = calculate_p(reactions, bacteria.rxn.matrix, 'Class', organism, scan_list, res_list)
    
  } else if (organism %in% bacteria.rxn.matrix$Phylum) {
    res = calculate_p(reactions, bacteria.rxn.matrix, 'Phylum', organism, scan_list, res_list)
    
  } else if (organism %in% bacteria.rxn.matrix$Domain) {
    res = no_optim(reactions, bacteria.rxn.matrix)
    
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Genus) {
    res = calculate_p(reactions, imgm.archaea.rxn.matrix, 'Genus', organism, scan_list, res_list)
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Family) {
    res = calculate_p(reactions, imgm.archaea.rxn.matrix, 'Family', organism, scan_list, res_list)
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Order) {
    res = calculate_p(reactions, imgm.archaea.rxn.matrix, 'Order', organism, scan_list, res_list)
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Class) {
    res = calculate_p(reactions, imgm.archaea.rxn.matrix, 'Class', organism, scan_list, res_list)
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Phylum) {
    res = calculate_p(reactions, imgm.archaea.rxn.matrix, 'Phylum', organism, scan_list, res_list)
    
  } else if (organism %in% imgm.archaea.rxn.matrix$Domain) {
    res = no_optim(reactions, imgm.archaea.rxn.matrix)
    
  } else {
    stop('Organism taxonomy not detected for ', organism, ', even at the domain level. Please see MetaPredict usage instructions [-h].') }
  
  return(res)
}



#function to reconstruct and predict KEGG metabolic pathways
metapredict = function(userData, orgNames) {
  orgs = read_csv(orgNames, col_names = 'col', col_types = cols())
  
  temp = userData %>%
    inner_join(ko_term.tibble, by = 'ko_term') %>%
    filter(!(duplicated(reaction))) #at this point we have ko_term, reaction, and reaction_description columns
  
  res = temp %>%
    group_split() %>%
    map(full_join, pathways.tibble, temp, by = c('reaction', 'reaction_description')) %>%
    map(arrange, pathway) %>%
    map(group_by, pathway) %>%
    map(filter, !(all(is.na(organism))), !(is.na(pathway)))
  
  scan_missing = res %>%
    map(ungroup) %>%
    map(filter, is.na(ko_term), reaction %in% colnames(bacteria.rxn.matrix)) %>% #bacteria/archaea matrices have the same column names
    map(select, -ko_term)
  
  map(1:length(orgs$col), ~ { #added scan here...changed 'predictions' to 'results'..
    #message('Processing ', orgs$col[.x], '...')
    message('Processing ', unique(userData$organism)[.x], '...')
    pull_data(scan_missing[[.x]]$reaction, orgs$col[[.x]], scan_missing[[.x]],
              res[[.x]]) }) # added .y (scan) here...
  #this can probably be combined with command below, using the ~ map syntax
  #pred = map(predictions, function(prediction) {
  #  tibble(reaction = as_vector(map(prediction, names)),
  #         probability = as_vector(prediction)) }) 
      #left_join(ko_term.tibble, by = 'reaction') %>%
      #select(-ko_term) })
  
  #scan = scan %>%
  #  map2(pred, ~ mutate(.x, probability = .y$probability)) %>%
  #  map(group_by, pathway)
  
  #out = map2(res, scan, full_join, by = c('reaction', 'reaction_description', 
  #                                        'pathway', 'pathway_name', 'pathway_class', 'organism')) %>%
  #  map(mutate, probability = as.character(probability)) %>%
  #  map(mutate, probability = case_when(!(is.na(ko_term)) & is.na(probability) ~ 'Present',
  #                              TRUE ~ probability))
  
  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')
  
  #return(results)  
}



#create options to parse user command line arguments
option_list = list(
  make_option(c('-p', '--path'), action = 'store',
              help = 'Path to directory containing KEGG orthology files', type = 'character',
              default = NA),
  
  make_option(c('-f', '--filePattern'), action = 'store',
              help = 'Regex suffix pattern for KEGG orthology files [default %default]',
              type = 'character', default = '*-ko.tsv'),
  
  make_option(c('-x', '--taxonList'),
              help = 'Path to a CSV file of genera of each organism, with each genus on a new line, 
              in the same order as the corresponding KEGG orthology files listed in --path',
              action = 'store', type = 'character', default = NA),
  
  make_option(c('-c','--cores'), default = 1, help = 'Number of CPU cores to utilize [default %default]', 
              action = 'store', type = 'integer'),
  
  make_option(c('-o', '--output'), action = 'store', default = './',
              help = 'Path to output directory. It will be created if it does not already exist [default %default]',
              type = 'character'),
  
  make_option(c('-v', '--verbose'), action = 'store_true', default = F,
              help = 'If TRUE, lists progress information for each reaction probability calculation, 
              prints one progress bar for entire job if FALSE [default %default]'),
  
  make_option(c('-e', '--e-value'), action = 'store', default = 1e-3,
              help = 'E-value for hmm hits from Kofamscan. All hits above the given value will be 
              discarded [default %default]', type = 'double')
)



#save user command line arguments to argv object
message('Parsing command line arguments...')
argv = parse_args(OptionParser(option_list = option_list))
message('Done\n')

if(all(!is.na(c(argv$path, argv$taxonList)))) {

#load required objects for analysis
message('Loading required data objects...')
load('../reqd-metapredict-data-objects-v3.RData')
message('Done\n')

#message('Connecting to SQL database...')
#meta.DB = dbConnect(duckdb::duckdb(), '/vortexfs1/omics/pachiadaki/dgellermcgrath/scripts/metapredict/sql_databases')
#message('Done\n')

#setting parallel computing settings based on user input
message('Setting parallel computing parameters...')
plan(multicore, workers = argv$cores)
options(future.globals.maxSize = 3145728000) #setting max mem for objects to 3GB [default is 500MB]
message('Done\n')

#read in user data, run MetaPredict
message('Parsing HMM hits and E-values into MetaPredict. Using E-value cutoff: ', argv$`e-value`)
userData = read_data(argv$path, argv$filePattern)
message('Done\n\nPerforming metabolic pathway reconstruction and reaction predictions...\n')
metapredict(userData, argv$taxonList)
#result = metapredict(userData, argv$genusList)


#save results to output folder, one TSV for each organism
#outpath = pmap(list(argv$output, unique(userData$organism), '-MetaPredict.tsv'), paste, sep = '')
#message('\nSaving output...')
#map2(result, outpath, write_tsv)
#dbDisconnect(meta.DB)
message('All done.')

} else {
    stop('You did not specify --path and --taxonList arguments properly. Please see usage')
}





#still under development: 
#function to parse in new KEGG Orthology data for genomes; still in development. this will allow users to
#add their own metabolic data from their own microbial genomes to MetaPredict, to include in its calculations
#add_new_data = function(filePath, filePattern, flatfile, flatDelim, workingDir) {
#  setwd(filePath)
#  index = list.files(path = filePath, pattern = filePattern)
#  res = list()

#NOTE: need to add in function to read in a req'd user flatfile which gives genus, a unique organism ID, and domain (bacteria or archaea) for each add
#  user_key = read_delim(flatfile, delim = flatDelim) %>%
#    mutate(colnames = case_when(
#      str_detect(colnames(user_key), regex('Genus', ignore_case = T)) ~ 'Genus',
#      str_detect(colnames(user_key), regex('Genome.ID', ignore_case = T)) ~ 'Genome.ID',
#      str_detect(colnames(user_key), regex('Domain', ignore_case = T)) ~ 'Domain',
#      TRUE ~ colnames(user_key))) %>%
#    select(Genus, Genome.ID, Domain)

#  future_map(1:length(index), .progress = T, ~ {
#    col = sub('(.*)-ko.tsv', '\\1', index[.x], perl = T)
#    res[[.x]] = read_delim(index[.x], col_types = cols(), delim = '\t') %>%
#      slice(-1) %>%
#      select(KO, `E-value`) %>%
#      mutate(`E-value` = as.numeric(`E-value`)) %>%
#      filter(`E-value` <= argv$`e-value`) %>% #e-value is set by argv, default = 1e-3
#      select(KO) %>%
#      rename(!!col := KO) %>%
#      filter(!(duplicated(.data[[col]])))
#  })
#  res = tibble(res)

#  res.rxn.matrix = future_map_dfc(1:length(res[[1]]), function(x)
#    map_lgl(1:length(rxn.tib[[1]]), function(y) 
#      any(str_detect(res[[1]][[x]][[1]], rxn.tib[[1]][[y]][[1]])))) %>%
#    rename_all(funs(map_chr(res[[1]], names))) %>%
#    map_dfc(as.numeric) %>%
#    add_column(col = map_chr(rxn.tib[[1]], names), .before = 1) %>%
#    column_to_rownames(var = 'col') %>%
#    t() %>%
#    as.data.frame() %>%
#    rownames_to_column(var = 'Genome.ID') %>%
#    inner_join(select(user_key, Genome.ID, Genus, Domain), by = 'Genome.ID') %>%
#    select(Genome.ID, Genus, Domain, everything()) %>%
#    arrange(Domain) %>%
#    group_by(Domain) %>%
#    group_split()

#  if (length(res.rxn.matrix == 2)) {  #need to verify that 'archaea' will always be grouped as the first of two lists
#    res.rxn.matrix = res.rxn.matrix %>% map(select, -Domain)

#    imgm.archaea.rxn.matrix = imgm.archaea.rxn.matrix %>%
#      ungroup() %>%
#      bind_rows(res.rxn.matrix[[1]]) %>%
#      group_by(Genus)

#    bacteria.rxn.matrix = bacteria.rxn.matrix %>%
#      ungroup() %>%
#      bind_rows(res.rxn.matrix[[2]]) %>%
#      group_by(Genus)

#  } else if (length(res.rxn.matrix == 1)) {
#    if (all(res.rxn.matrix[[1]]$Domain == regex('bacteria', ignore_case = T))) {
#      res.rxn.matrix = res.rxn.matrix %>% map(select, -Domain)

#      bacteria.rxn.matrix = bacteria.rxn.matrix %>%
#        ungroup() %>%
#        bind_rows(res.rxn.matrix[[1]]) %>%
#        group_by(Genus)

#    } else if (all(res.rxn.matrix[[1]]$Domain == regex('archaea', ignore_case = T))) {
#      res.rxn.matrix = res.rxn.matrix %>% map(select, -Domain)

#      imgm.archaea.rxn.matrix = imgm.archaea.rxn.matrix %>%
#        ungroup() %>%
#        bind_rows(res.rxn.matrix[[1]]) %>%
#        group_by(Genus)

#    } else {
#      stop("The domain is not specified for one or more genomes as either 'Archaea' or 'Bacteria' (case insensitive).") }

#  } else {
#    stop('Please check that your input flatfile contains Domain, Genus, and Genome.ID columns') }

#  save(bacteria.rxn.matrix, imgm.archaea.rxn.matrix, pathways.tibble, ko_term.tibble,
#       file = paste(workingDir, 'reqd-metapredict-data-objects.RData'))
#}

