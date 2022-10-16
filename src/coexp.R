#Co-expression analyses - using reference co-expression networks from CoCoCoNet (Lee et al., 2020, NAR)

library(rhdf5)

#Load in list of species and taxonomical distance matrix across species

list_species <- c("human", "rat", "mouse", "cow", "boar", "bee", "rhesusm", "sheep", "goat", "chicken", "crabm", "dog", "horse", "zebrafish", "chimp")
species <- read.csv("species_classification_ranks_processed.txt")
taxo_distances <- read.csv("spe39_divergence_timetree.csv")

#Load in ortholog maps to find 1-1 orthologs of cluster-specific disease signatures across chordate species (from CoCoCoNet - Lee et al., 2020)

all_maps <- list.files(path = "CoCoCoNet/orthologMaps/")
all_maps <- data.frame(all_maps) %>%
  mutate(first = str_split(all_maps, "_", simplify = TRUE)[, 1], second = str_split(all_maps, "_", simplify = TRUE)[, 2]) %>%
  filter(first == "human" & second %in% list_species | first %in% list_species & second == "human")

#Load in all reference co-expression networks except the human network and convert row and column names

chordata <- species %>% filter(saved_name != "human") %>% filter(phylum == "Chordata") %>% dplyr::select(saved_name) %>% "[["(1)

all_maps <- list.files(path = "CoCoCoNet/orthologMaps/")
all_maps <- data.frame(all_maps) %>%
  mutate(first = str_split(all_maps, "_", simplify = TRUE)[, 1], second = str_split(all_maps, "_", simplify = TRUE)[, 2]) %>%
  filter(first == "human" & second %in% chordata | first %in% chordata & second == "human")

ref_networks <- list()

z <- lapply(chordata, FUN = function(i){
  print(i)
  species_net_path <- paste("CoCoCoNet/networks/", i, "_prioAggNet.hdf5", sep = "")
  col <- h5read(species_net_path, "col")
  row <- h5read(species_net_path, "row")
  species_net <- h5read(species_net_path, "agg")
  rownames(species_net) <- row
  colnames(species_net) <- col
  ref_networks[[i]] <<- species_net
})

#Load in human co-expression network which is saved as a HDF5 file

cococonet_matrix <- "CoCoCoNet/networks/human_prioAggNet.hdf5"
col <- h5read(cococonet_matrix, "col")
row <- h5read(cococonet_matrix, "row")
dat <- h5read(cococonet_matrix, "agg")
rownames(dat) <- str_split(row, "\\.", simplify = TRUE)[, 1]
colnames(dat) <- str_split(col, "\\.", simplify = TRUE)[, 1]

#Add human reference co-expression network to list of other networks and free up memory

ref_networks[['human']] <- dat
rm(dat)
