library('tidyverse')
library('visNetwork')
library('htmlwidgets')
library('igraph')
library('Biostrings')
library('ggtree')
library('ape')

system.file("tex", "texshade.sty", package="msa")

local <- Biostrings::readDNAStringSet('local_db_solo/multifasta_x3.fasta')
ncbi <- Biostrings::readDNAStringSet('ncbi_repeats_db/ncbi_repeats_x3.fasta')
both <- c(local, ncbi)

blastres <- read_csv('blastres.csv')

bestblastres <- blastres %>% 
  dplyr::filter(sum_coverage >= 0.7) %>%
  dplyr::group_by(qseqid, sseqid) %>% 
  dplyr::arrange(desc('sum_coverage'), 'task') %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(stringr::str_detect(db, '^loca.*')) %>% 
  dplyr::select(qseqid, sseqid) %>% 
  dplyr::rename(from=qseqid, to=sseqid)

blast_graph <- bestblastres %>% 
igraph::graph_from_data_frame()

dg <- igraph::decompose.graph(blast_graph)

names_list <- dg %>% map(~{.x %>% V %>% .$name})

names_list %>% 
  purrr::map(~ paste(.x, collapse='\t')) %>% 
  unlist() %>% 
  paste0(collapse = '\n') %>% 
  write_file('clusters/list_of_clusters.tsv')

slice_fasta_by_name_prefix <- function(x, both){
  `%>%` <- magrittr::`%>%`
  map(x, ~both@ranges %>% as.data.frame() %>% tibble::as_tibble() %>% 
    tibble::rownames_to_column('id') %>% 
    dplyr::filter(str_starts(names, .x)) %>%
    dplyr::mutate(id = as.integer(id)) %>%
    dplyr::pull(id)) %>%
    unlist()
}

all_groups <- names_list %>% 
  purrr::map(\(x){both[slice_fasta_by_name_prefix(x, both)]})

all_groups %>% 
  purrr::walk2(., seq_along(.),
    ~{Biostrings::writeXStringSet(.x, 
                                glue::glue("clusters/{.y}.fasta"))})


Sys.glob('clusters/*.fasta.maf') %>% purrr::walk(~{

tryCatch({
  nbin <- Biostrings::readDNAMultipleAlignment(.x) %>%
    ape::as.DNAbin() %>%
    ips::trimEnds(min.n.seq = length(.))
  
  old_names <- nbin %>% 
    ape::as.alignment() %>% 
    {.$nam}
  
  num_bp = old_names %>% 
    stringi::stri_match_first_regex('\\d+(?=nt)') %>% 
    as.numeric() %>% 
    {max(.) * 1.1} %>% 
    floor()
  
  new_names <- old_names %>% 
    purrr::map_vec(~stringr::str_replace_all(.x, '(_R_)|(_TR_.*)', ''))
  
  nbin <- nbin %>% 
    ape::as.alignment()
  
  nbin$nam <- new_names
  
  nbin <- nbin %>% 
    ape::as.DNAbin() %>% 
    {.[,1:num_bp]}
  
  tree <- nbin %>%
    ape::dist.dna() %>%
    ape::nj() %>% 
    ape::multi2di()
  
  tree$edge.length <- pmax(tree$edge.length,0.001)
  
  boots <- ape::boot.phylo(tree, nbin, function(xx) nj(dist.dna(xx)), 
                      B = 10000,  mc.cores = 6) %>% 
    {round(100*./10000)}
  
  bs_tibble_tip <- tibble(
    node=1:Nnode(tree) + Ntip(tree), 
    bootstrap = ifelse(boots < 50, "", boots))
  
  bs_tibble_node <- tibble(
    node=1:Ntip(tree), 
    source = tree$tip.label) %>% 
    mutate(source =  stringi::stri_replace_first_regex(source, '_.*', '')) %>% 
    mutate(color = case_when(
      source %in% c('V5', 'KP1', 'KP2') ~ 'P. spicata',
      source %in% c('KP3', 'KP4', 'KP5') ~ 'P. libanotica',
      source %in% c('KP10') ~ 'P. tauri',
      source %in% c('KA25', 'KA26', 'V1') ~ 'Th. bessarabicum',
      source %in% c('KK3') ~ 'Ae. tauschii typica',
      source %in% c('KK1') ~ 'Ae. tauschii strangulata',
      source %in% c('KK2') ~ 'Ae. tauschii',
      source %in% c('KK6') ~ 'Ae. crassa 4x',
      source %in% c('KA27') ~ 'Ae. crassa 6x'
    )) %>% 
    mutate(name = glue::glue("{source} {color}"))
  
  bs_tibble <- full_join(bs_tibble_tip, bs_tibble_node)
  
  
  tree_plot <- tree %>%
    ggtree(branch.length = 'none') %<+% bs_tibble +
    #geom_text(, hjust=-.25, size = 3) +
    geom_nodelab(aes(label = bootstrap), size = 5, nudge_x = 0.25) +
    geom_tiplab(aes(label = name, color = color), align=TRUE, size = 5) +
    geom_tippoint(aes(color = color))


msaplot(tree_plot, nbin, offset=2.5, width=1, height = 0.5) -> plot_to_save

ggsave(glue::glue("{.x}.png"), plot=plot_to_save, width = 20, height = 12)
},
error=print)
}
)
  # as.matrix() %>%
  # ips::trimEnds() %>%
  # ape::as.alignment()

saveWidget(visIgraph(dg[[1]]), file = "cluster_interactive_example.html")



Sys.glob('rep_*.xlsx')

Sys.glob('rep_*.xlsx') %>% 
  purrr::map_dfr(
    ~ readxl::read_excel(.x, sheet = 'Sheet1') %>% 
      janitor::clean_names() %>% 
      select(x1, sum_percent) %>% 
      mutate(source = .x  %>% 
               stringi::stri_extract_first_regex('[A-Z]+\\d+'))
  ) %>% 
  pivot_wider(names_from = x1, values_from = sum_percent) %>% 
  janitor::remove_constant() %>% 
  mutate(across(
    where(is.numeric), scale
  )) -> structure_table

dist_matrix <- proxy::dist(structure_table %>% select(-source), 
                           diag = TRUE, pairwise = TRUE) %>% 
  as.matrix()

row.names(dist_matrix) <- c(structure_table %>% pull(source))
colnames(dist_matrix) <- c(structure_table %>% pull(source))

tree <- dist_matrix %>% ape::njs() %>% ape::multi2di()
  
tree$edge.length <- pmax(tree$edge.length,0.001)

boots <- ape::boot.phylo(tree, dist_matrix, function(xx) njs(xx), 
                         B = 10000,  mc.cores = 6) %>% 
  {round(100*./10000)}

bs_tibble_tip <- tibble(
  node=1:Nnode(tree) + Ntip(tree), 
  bootstrap = ifelse(boots < 50, "", boots), boots))

bs_tibble_node <- tibble(
  node=1:Ntip(tree), 
  source = tree$tip.label) %>% 
  mutate(source =  stringi::stri_replace_first_regex(source, '_.*', '')) %>% 
  mutate(color = case_when(
    source %in% c('V5', 'KP1', 'KP2') ~ 'P. spicata',
    source %in% c('KP3', 'KP4', 'KP5') ~ 'P. libanotica',
    source %in% c('KP10') ~ 'P. tauri',
    source %in% c('KA25', 'KA26', 'V1') ~ 'Th. bessarabicum',
    source %in% c('KK3') ~ 'Ae. tauschii typica',
    source %in% c('KK1') ~ 'Ae. tauschii strangulata',
    source %in% c('KK2') ~ 'Ae. tauschii',
    source %in% c('KK6') ~ 'Ae. crassa 4x',
    source %in% c('KA27') ~ 'Ae. crassa 6x'
  )) %>% 
  mutate(name = glue::glue("{source} {color}"))

bs_tibble <- full_join(bs_tibble_tip, bs_tibble_node)

tree_plot <- tree %>%
  ggtree(branch.length = 'none') %<+% bs_tibble +
  geom_nodelab(aes(label = bootstrap), size = 5, nudge_x = 0.25) +
  geom_tiplab(aes(label = name, color = color), align=TRUE, size = 5, offset = 0.1) +
  geom_tippoint(aes(color = color)) + hexpand(1)

tree_plot

ggsave(glue::glue("structure.png"), plot=tree_plot, width = 20, height = 12)

library(ggridges)

qpcr_csv <- read_tsv('qpcr.tsv') %>% 
  pivot_longer(cols = contains('.'), 
               names_to = 'Вид', 
               values_to = 'Log10 от доли относительно VRN1') %>% 
  mutate(cl_num = stringi::stri_extract_first_regex(Повтор, '\\d+') %>% as.numeric) %>% 
  arrange(cl_num, Вид)

ggplot(qpcr_csv, mapping = aes(x = Повтор, y = `Log10 от доли относительно VRN1`,
                                 fill = Вид)) +
  geom_bar(stat = "identity", position = "dodge")


