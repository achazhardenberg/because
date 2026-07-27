library(because)
library(because.phybase)
library(ape)
library(MASS)

set.seed(42)
BETA_BM_MR          <-  0.80   
BETA_MR_TT          <-  0.00   
BETA_BM_TT          <-  0.60   
BETA_ELEV_URES      <-  0.00   
BETA_URES_NDVI      <-  0.50   
BETA_ELEV_NDVI      <-  0.00   
BETA_URES_FC        <-  0.60   
BETA_ELEV_FC        <-  0.20   
BETA_ELEV_TEMP      <- -1.50   
SD_TEMP_WITHIN      <-  1.00   
BETA_FC_ABUND       <-  0.50   
BETA_WIND_ABUND     <- -0.20   
BETA_MISMATCH_ABUND <- -0.50   
PAGEL_LAMBDA        <-  0.60   

N_Species          <- 50
N_Site             <- 30   
N_Surveys_per_Site <-  3   
N_Surveys          <- N_Site * N_Surveys_per_Site   

species_tree <- rtree(N_Species)
species_tree$tip.label  <- paste0("Sp_", 1:N_Species)
species_tree$edge.length <- species_tree$edge.length /
    max(node.depth.edgelength(species_tree))

Vcv      <- vcv.phylo(species_tree)
V_lambda <- PAGEL_LAMBDA * Vcv + (1 - PAGEL_LAMBDA) * diag(diag(Vcv))

Body_Mass_s <- as.vector(mvrnorm(1, rep(0, N_Species), V_lambda))
eps_MR      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.15^2 * V_lambda))
Metabolic_Rate <- BETA_BM_MR * Body_Mass_s + eps_MR
eps_TT      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.40^2 * V_lambda))
Thermal_Tol <- BETA_BM_TT * Body_Mass_s + eps_TT 
names(Body_Mass_s) <- names(Metabolic_Rate) <- names(Thermal_Tol) <-
    species_tree$tip.label

d_species <- data.frame(
    Species        = species_tree$tip.label,
    Body_Mass_s    = Body_Mass_s,
    Metabolic_Rate = Metabolic_Rate,
    Thermal_Tol    = Thermal_Tol
)

site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
V_spatial <- exp(-as.matrix(dist(site_coords)) / 1.0)

Elevation_s <- as.vector(MASS::mvrnorm(1, rep(0, N_Site), V_spatial))
U_Resource   <- rnorm(N_Site, 0, 1.0)
NDVI         <- BETA_URES_NDVI * U_Resource +
                BETA_ELEV_NDVI * Elevation_s +
                rnorm(N_Site, 0, 0.20)
Flower_Cover <- BETA_URES_FC   * U_Resource +
                BETA_ELEV_FC   * Elevation_s +
                rnorm(N_Site, 0, 0.20)
d_site <- data.frame(
    Site         = paste0("Site_", 1:N_Site),
    Elevation_s  = Elevation_s,
    NDVI         = NDVI,
    Flower_Cover = Flower_Cover
)

site_of_survey  <- rep(1:N_Site, each = N_Surveys_per_Site)
Temperature <- BETA_ELEV_TEMP * Elevation_s[site_of_survey] +
               rnorm(N_Surveys, 0, SD_TEMP_WITHIN)
Wind_Speed  <- rnorm(N_Surveys, 0, 1)
d_survey <- data.frame(
    Survey      = paste0("Survey_", 1:N_Surveys),
    Site        = d_site$Site[site_of_survey],
    Temperature = Temperature,   
    Wind_Speed  = Wind_Speed     
)

d_obs <- expand.grid(
    Survey  = d_survey$Survey,
    Species = d_species$Species,
    stringsAsFactors = FALSE
)
d_obs <- merge(d_obs, d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed")],
               by = "Survey")
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs <- merge(d_obs, d_site[, c("Site", "Flower_Cover")], by = "Site")
d_obs$obs_id <- seq_len(nrow(d_obs))

Thermal_Mismatch <- d_obs$Temperature - d_obs$Thermal_Tol
u_site    <- rnorm(N_Site, 0, 0.40)
u_survey  <- rnorm(N_Surveys, 0, 0.40)
u_species <- rnorm(N_Species, 0, 0.40)
names(u_site)    <- paste0("Site_", 1:N_Site)
names(u_survey)  <- paste0("Survey_", 1:N_Surveys)
names(u_species) <- paste0("Sp_", 1:N_Species)

log_lambda <- BETA_FC_ABUND       * d_obs$Flower_Cover   +
              BETA_WIND_ABUND     * d_obs$Wind_Speed      +
              BETA_MISMATCH_ABUND * Thermal_Mismatch      +
              u_site[d_obs$Site]                          +
              u_survey[d_obs$Survey]                      +
              u_species[d_obs$Species]                    +
              rnorm(nrow(d_obs), 0, 0.10)
d_obs$Abundance <- rpois(nrow(d_obs), exp(log_lambda))

data_list <- list(
    species = d_species,
    site    = d_site,
    survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed")],
    obs     = d_obs[, c("obs_id", "Survey", "Species", "Abundance")]
)

structures <- list(
    phylo   = species_tree,   
    spatial = V_spatial        
)

eqs <- list(
    NDVI         ~ U_Resource,
    Flower_Cover ~ U_Resource + Elevation_s,
    Body_Mass_s    ~ 1,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol    ~ Body_Mass_s,
    Temperature  ~ Elevation_s,
    Abundance ~ Flower_Cover + Wind_Speed +
                I(Temperature - Thermal_Tol) +
                (1 | Site) + (1 | Survey)
)

print("Starting because...")
fit_dsep <- because(
    eqs,
     data      = data_list,
     levels    = list(
         species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
         site    = c("Elevation_s", "NDVI", "Flower_Cover", "U_Resource", "Site"),
         survey  = c("Temperature", "Wind_Speed", "Survey"),
         obs     = c("Abundance")),
     hierarchy    = "site > survey > obs; species > obs",
     link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
     latent       = "U_Resource",
     latent_method = "explicit",
     family       = c(Abundance = "negbinomial"),
     structure    = structures,
     dsep         = TRUE,     
     parallel     = TRUE,
     n.cores = 3,
     engine = "numpyro")
