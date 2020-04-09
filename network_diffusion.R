# network_diffusion.R
# uses STRINGdb to get a human PPI as an igraph object, adds some edges to
# COVID proteins, then does graph diffusion to see what other human
# proteins are "close"

library(STRINGdb)
library(igraph)
library(tidyverse)

string_db <- STRINGdb$new(version = "10", species = 9606)
string_net <- string_db$get_graph()
id_lookup <- read.csv(
  "string_ids_for_covid_interactors.csv",
  header = F,
  col.names = c("uniprot", "string")  
)
interactors_raw <- read.csv("supp_table_2.csv", header = F, sep = ",", skip = 1)
interactors <- interactors_raw %>%
  select(V1, V8)
interactors
