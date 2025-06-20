data {
  int<lower=1> num_indiv;
  int<lower=1> num_types;
  int<lower=1> num_pot;
  int<lower=0> num_pt_groups;
  int<lower=0> n_cells;
  int<lower=0> d_cells;
  int<lower=0> grid_size; // number of grid points
  vector[n_cells] oset;


  int<lower=1> n_samples;
  array[n_samples,2] int<lower=1,upper=n_cells> y_start_stop;

  matrix[n_cells,d_cells] x_cells;
}

parameters {
  matrix[d_cells,num_pt_groups] beta_global;
  matrix[d_cells,num_indiv] beta_indiv;
  array[n_samples] vector[d_cells] beta_local;
  vector<lower=0>[num_pt_groups > 0 ? num_pot : 0] sigma_beta_global;
  vector<lower=0>[num_pt_groups > 0 ? 1 : num_pot] sigma_beta_indiv;
  real<lower=0> sigma_beta_local;
  vector<lower=0>[num_pt_groups > 0 ? 1 : 0] tau_alpha_global;
  real<lower=0> tau_alpha_indiv;
  real<lower=0> tau_alpha_local;
}

generated quantities {
   array[n_samples] vector[grid_size] p_pred;
   {
    for (s in 1:n_samples) {
      int y_start = y_start_stop[s,1];
      int y_stop = y_start_stop[s,2];

      p_pred[s] = x_cells[y_start:y_stop] * beta_local[s] + oset[y_start:y_stop];

    }
  }


}
