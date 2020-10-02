#load required packages and their dependencies
if(all(c('tidyverse', 'furrr', 'optparse', 'progress') %in% rownames(installed.packages()))) {
  invisible(lapply(c('tidyverse', 'furrr', 'optparse', 'progress'), function(package) {
    suppressWarnings(suppressPackageStartupMessages(library(package, character.only = T)))})) 
  } else {
  stop('One or more of the required packages are not installed. Are you currently running the 
       MetaPredict conda environment?')
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



#function to parse in user data for pathway reconstruction & reaction prediction
read_data = function(filePath, filePattern) {
  setwd(filePath)
  index = list.files(path = filePath, pattern = filePattern)
  res = list()
  
  future_map(1:length(index), .progress = T, ~ {
    col = sub('(.*)-ko.tsv', '\\1', index[.x], perl = T)
    res[[.x]] = read_delim(index[.x], col_types = cols(), delim = '\t') %>%
      slice(-1) %>%
      select(KO, `E-value`) %>%
      mutate(`E-value` = as.numeric(`E-value`)) %>%
      filter(`E-value` <= argv$`e-value`) %>% #e-value is set by argv, default = 1e-3
      select(KO) %>%
      rename(!!col := KO) %>%
      filter(!(duplicated(.data[[col]])))
  })
  
  res = res %>%
    tibble() %>%
    unnest(cols = everything()) %>% 
    #gather(key = 'organism', value = 'ko_term', na.rm = T) %>%
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



#function to find maximum likelihood estimates of alpha and beta
calculate = function(reactions, organism)  {
  #if (argv$verbose == T) {
  #  pb = progress_bar$new(show_after = 0, total = length(reactions), clear = F,
  #                        format = '[:bar] :current/:total | :percent complete | eta: :eta | time elapsed: :elapsedfull',
  #                        callback = invisible, width = 200)
  #} else {
  #  pb = progress_bar$new(show_after = 0, total = length(reactions), clear = F,
  #                        format = 'Calculating :what [:bar] :current/:total | :percent complete | eta: :eta | time elapsed: :elapsedfull',
  #                        callback = invisible, width = 200) }
  
  future_map(reactions, .progress = T, ~ {         #need to play with the future_map built-in progress bar...
    if (organism %in% bacteria.rxn.matrix$Genus) {
      coll.k = bacteria.rxn.matrix %>%
        filter(Genus != organism) %>%
        summarize(y_k = sum(.data[[.x]]), n_k = length(Genus)) 
      
      coll.j = bacteria.rxn.matrix %>%
        filter(Genus == organism) %>%
        summarize(y_j = sum(.data[[.x]]), n_j = length(Genus))
      
    } else if (organism %in% imgm.archaea.rxn.matrix$Genus) {
      coll.k = imgm.archaea.rxn.matrix %>%
        filter(Genus != organism) %>%
        summarize(y_k = sum(.data[[.x]]), n_k = length(Genus)) 
      
      coll.j = imgm.archaea.rxn.matrix %>%
        filter(Genus == organism) %>%
        summarize(y_j = sum(.data[[.x]]), n_j = length(Genus))
      
    } else {
      stop('Genus not found. Please see MetaPredict usage instructions [-h].')
    }

   # if (argv$verbose == T) {
  #  message('Calculating alpha and beta maximum likelihood estimates for reaction ', 
  #        .x, ' for organism ', organism, '...', sep = '') 
  #    pb$tick()
  #  
  #    } else {
  #      pb$tick(tokens = list(what = paste(organism, ': ', .x, sep = ''))) }
    
    #opt = optim(par = c(0.01, 0.01), LogL.bb.5, y_k = coll.k$y_k, n_k = coll.k$n_k)
                #method = 'L-BFGS-B', lower = 0.01, upper = 10000)
    
    opt = optim(par = c(0.01, 0.01), LogL.bb.5, y_k = coll.k$y_k, n_k = coll.k$n_k,
                method = 'L-BFGS-B', lower = 1e-10, upper = 1000000)
    
    res = (opt$par[1] + coll.j$y_j) / (opt$par[1] + opt$par[2] + coll.j$n_j)
    names(res) = .x
    
    if (argv$verbose == T) {
    message('\nDone\n') }
    
    return(res) })
}



#function to reconstruct and predict KEGG metabolic pathways
metapredict = function(userData, orgNames) {
  orgs = read_csv(orgNames, col_names = 'col', col_types = cols())
  
  temp = userData %>% 
    inner_join(ko_term.tibble, by = 'ko_term') %>%
    filter(!(duplicated(reaction)))
  
  res = temp %>%
    group_split() %>%
    map(full_join, pathways.tibble, temp, by = 'reaction') %>%
    map(arrange, pathway) %>%
    map(group_by, pathway) %>%
    map(filter, !(all(is.na(organism))), !(is.na(pathway)))
  
  scan = res %>%
    map(ungroup) %>%
    map(filter, is.na(ko_term), reaction %in% colnames(bacteria.rxn.matrix)) %>% #bacteria/archaea matrices have the same column names
    map(select, reaction, pathway)
  
  predictions = map(1:length(orgs$col), ~ {
    calculate(scan[[.x]]$reaction, orgs$col[[.x]]) })
  
  #this can probably be combined with command below, using the ~ map syntax
  pred = map(predictions, function(prediction) {
    tibble(reaction = as_vector(map(prediction, names)),
           probability = as_vector(prediction)) })
  
  scan = scan %>% 
    map2(pred, ~ mutate(.x, probability = .y$probability)) %>%
    map(group_by, pathway)
  
  out = map2(res, scan, full_join, by = c('reaction', 'pathway'))
  
  message('Finished KEGG metabolic pathway reconstruction and reaction probability calculations.')
  
  return(out)  
}



#create options to parse user command line arguments
option_list = list(
  make_option(c('-p', '--path'), action = 'store',
              help = 'Path to directory containing KEGG orthology files', type = 'character',
              default = NA),
  
  make_option(c('-f', '--filePattern'), action = 'store',
              help = 'Regex suffix pattern for KEGG orthology files [default %default]',
              type = 'character', default = '*-ko.tsv'),
  
  make_option(c('-g', '--genusList'),
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

if(all(!is.na(c(argv$path, argv$genusList)))) {

#load required objects for analysis
message('Loading required data objects...')
load('../reqd-metapredict-data-objects.RData')
message('Done\n')

#setting parallel computing settings based on user input
message('Setting parallel computing parameters...')
plan(multicore, workers = argv$cores)
options(future.globals.maxSize = 3145728000) #setting max mem for objects to 3GB [default is 500MB]
message('Done\n')

#read in user data, run MetaPredict
message('Parsing HMM hits and E-values into MetaPredict. Using E-value cutoff: ', argv$`e-value`)
userData = read_data(argv$path, argv$filePattern)
message('Done\n\nPerforming metabolic pathway reconstruction and reaction predictions...')
result = metapredict(userData, argv$genusList)

#save results to output folder, one TSV for each organism
outpath = pmap(list(argv$output, unique(userData$organism), '-MetaPredict.tsv'), paste, sep = '')
map2(result, outpath, write_tsv)
message('\nAll done. Saving output.')

} else {
    stop('You did not specify --path and --genusList arguments properly. Please see usage')
}
