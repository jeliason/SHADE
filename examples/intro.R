library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)

devtools::load_all()

df 

df_raw %>%
  count(Spot,type) %>%
  group_by(type) %>%
  summarise(mn = median(n)) %>%
  filter(mn > 30) %>%
  pull(type) %>%
  as.character() -> types_keep

spots_pt <- pt_data

spots <- spots_pt %>%
  pull(Spot)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  path <- "PATH"
  # types_keep <- types_keep[1:2]
  # scale_factor <- 1
  n_dummy <- 1e3
  spots <- spots[c(1,2,5,6)]
  
} else {
  path <- "PATH"
  # scale_factor <- 2
  n_dummy <- 1e3
  args <- commandArgs(trailingOnly=TRUE)
  type_idx <- as.numeric(args[1])
}

spots_pt <- spots_pt %>%
  filter(Spot %in% spots)

num_types <- length(types_keep)

# file_json <- paste0(path,"data_stan_CRC_scale_factor_",scale_factor,".json")
# file_dlist <- paste0(path,"data_lists_CRC_scale_factor_",scale_factor,".rds")

num_pot <- 3
potentials <- make_rbfs(n_basis_functions = num_pot, max_dist = 75, basis_function_sigma = 15)
grainsize <- 1
mean_alpha <- -10
scale_sigmas <- 5

sample_to_indiv <- spots_pt %>%
  pull(Patient) %>%
  as.factor() %>%
  as.numeric()
indiv_to_group <- spots_pt %>%
  distinct(Patient,Group) %>%
  pull(Group) %>%
  as.factor() %>%
  as.numeric()

dats <- df_raw %>%
  dplyr::filter(Spot %in% spots) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels()
  })

pats <- lapply(1:length(dats),\(i) {
  df <- dats[[i]]
  tryCatch({
    pat <- make_pat(df$X,df$Y,factor(df$type,levels=types_keep))
    sq_W <- owin(xrange = c(min(df$X),max(df$X)),yrange=c(min(df$Y),max(df$Y)))
    Window(pat) <- sq_W
    pat
  },error=\(e) {
    print(i)
    stop(e)
  })
})

Qs <- lapply(pats,\(pat) make_quadrature(pat,n_dummy = n_dummy))

type <- types_keep[type_idx]
file_json <- paste0(path,"data_stan_CRC_type_",make.names(type),".json")
file_dlist <- paste0(path,"data_lists_CRC_type_",make.names(type),".rds")
file_metadata <- paste0(path,"metadata_CRC_type_",make.names(type),".rds")

lapply(Qs,\(Q) {
  log(intensity(Q$dummy)) %>%
    enframe() %>%
    right_join(tibble(name=marks(Q)),by="name") %>%
    filter(name == type) %>%
    pull(value)
}) %>% unlist() -> offset

data_lists <- lapply(1:length(Qs),\(i) {
  data_list <- make_acyclic_data(Qs[[i]],potentials,type,verbose=FALSE)
  print(i)
  data_list
})

saveRDS(data_lists,file_dlist)

x_cells <- do.call(rbind,lapply(data_lists,\(data_list) data_list$data))
# x_cells <- x_cells[,-1]

is_cell <- unlist(lapply(data_lists,\(data_list) data_list$response))
n_cells = nrow(x_cells)
d_cells = ncol(x_cells)
n_samples <- length(data_lists)
sample_id <- lapply(1:length(data_lists),\(i) {
  rep(i,length(data_lists[[i]]$response))
}) %>% unlist()

start_stop_idx <- lapply(1:length(data_lists),\(i) {
  idx <- sort(which(sample_id == i))
  c(idx[1],idx[length(idx)])
}) %>% do.call(rbind,.)

x_cells <- as(x_cells,"RsparseMatrix")

u <- as.integer(x_cells@p + 1)
w <- x_cells@x
v <- x_cells@j + 1

y_start_stop <- start_stop_idx

data_start_stop <- lapply(1:n_samples,\(i) {
  row_start <- y_start_stop[i,1]
  row_stop <- y_start_stop[i,2]
  
  data_start <- u[row_start]
  data_end <- u[row_stop + 1] - 1
  c(data_start,data_end)
}) %>% do.call(rbind,.)

num_indiv <- length(unique(spots_pt$Patient))
num_pt_groups <- length(unique(spots_pt$Group))

data_stan <- list(
  num_indiv = num_indiv,
  num_types = num_types,
  num_pot = num_pot,
  num_pt_groups = num_pt_groups,
  n_cells = n_cells,
  d_cells = d_cells,
  is_cell = is_cell,
  oset = -offset,
  n_samples = n_samples,
  sample_to_indiv = sample_to_indiv,
  indiv_to_group = indiv_to_group,
  sample_id = sample_id,
  
  y_start_stop = y_start_stop,
  data_start_stop = data_start_stop,
  n_nz = length(w),
  w = w,
  v = v,
  u = u,
  
  mean_alpha = mean_alpha,
  scale_sigmas = scale_sigmas,
  scale_sigma_betas = seq(5,1,length.out=length(potentials)),
  scale_sigma_alpha = scale_sigmas,
  
  grainsize = grainsize
)

print("writing json...")
write_json_chunked(data_stan, file_json,chunk_size = 1e6)
print("json written.")

metadata <- list(
  coef_names = colnames(x_cells),
  potentials = potentials,
  spots = unique(pt_data$Spot),
  types = types_keep,
  sample_to_indiv = sample_to_indiv,
  indiv_to_group = indiv_to_group
)

saveRDS(metadata,file_metadata)