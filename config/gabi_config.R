
######################################
###            METADATA            ###
######################################
# gen3sis configuration
#
# Version: 1.0
#
# Author:
#
# Date: 14.06.2024
#
# Landscape:
#
# Publications:
#
# Description: 
#
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation.
random_seed <- 20240614

# set the starting time step or leave NA to use the earliest/highest time-step.
start_time <- NA

# set the end time step or leave as NA to use the latest/lowest time-step (0).
end_time <- NA

# maximum total number of species in the simulation before it is aborted.
max_number_of_species <- 25000

# maximum number of species within one cell before the simulation is aborted.
max_number_of_coexisting_species <- 2500

# a list of traits to include with each species
# a "dispersal" trait is implicitly added in any case
trait_names <- c("dispersal", "biome", "continent")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges <- list()


######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation and collect additional information if desired.
end_of_timestep_observer <- function(data, vars, config) {
  plot_richness(data$all_species, data$landscape)
  # the list of all species can be found in data$all_species
  # the current landscape can be found in data$landscape
  # saving functions example:
    # save_landscape()
    # save_species()
  
  # plotting functions example:
    # plot environmental conditions
    # plot_landscape(data$landscape)
    # plot richness
    # plot_richness(data$all_species, data$landscape)
    # plot a specific environmental condition
    # plot_raster_single(data$landscape$environment[,"temp"], data$landscape, "temp", NA)
    # plot species 1 range
    # plot_species_presence(data$all_species[[1]], data$landscape)
    # plot(0,type="n",axes=FALSE,ann=FALSE)

}


######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when 
# colonizing a cell during the dispersal.
initial_abundance <- 1

# place species in the landscape:
create_ancestor_species <- function(landscape, config) {
  co <- landscape$coordinates
  biome <- landscape$environment[, "biome"]
  # Arid biomes
  starting_biome <- 2

  # One species starting in the south, one in the north
  sel_south <- co[, "x"] >= -70 &
    co[, "x"] <= -67 &
    co[, "y"] >= -42 &
    co[, "y"] <= -40 &
    biome == starting_biome
  cell_ids_south <- rownames(co)[sel_south]
  sp_south <- create_species(cell_ids_south, config)
  sp_south$traits[, "dispersal"] <- 1
  sp_south$traits[, "biome"] <- starting_biome
  sp_south$traits[, "continent"] <- 2

  sel_north <- co[, "x"] >= -104 &
    co[, "x"] <= -102 &
    co[, "y"] >= 31 &
    co[, "y"] <= 32 &
    biome == starting_biome
  cell_ids_north <- rownames(co)[sel_north]
  sp_north <- create_species(cell_ids_north, config)
  sp_north$traits[, "dispersal"] <- 1
  sp_north$traits[, "biome"] <- starting_biome
  sp_south$traits[, "continent"] <- 1

  list(sp_south, sp_north)
}


######################################
###             Dispersal          ###
######################################

# the maximum range to consider when calculating the distances from local distance inputs.
max_dispersal <- Inf

# returns n dispersal values.
get_dispersal_values <- function(n, species, landscape, config) {
  rweibull(n, shape = 1.5, scale = 200)
}


######################################
###          Speciation            ###
######################################

# threshold for genetic distance after which a speciation event takes place.
divergence_threshold <- 500

# factor by which the divergence is increased between geographically isolated population.
# can also be a matrix between the different population clusters.
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  1
}


######################################
###            Evolution           ###
######################################

# mutate the traits of a species and return the new traits matrix.
apply_evolution <- function(species, cluster_indices, landscape, config) {
  species$traits
}


######################################
###             Ecology            ###
######################################

# called for every cell with all occurring species, this function calculates abundances and/or 
# who survives for each sites.
# returns a vector of abundances.
# set the abundance to 0 for every species supposed to die.
apply_ecology <- function(abundance, traits, environment, config) {
  abundance[traits[, "biome"] != environment[, "biome"]] <- 0
  abundance
}
