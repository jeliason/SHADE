functions {
	real partial_sum_lpmf(array[] int slice_samples,
                        int sample_start_no_use, int sample_end_no_use,
												array[] int is_cell,
												array[,] int y_start_stop,
												array[,] int data_start_stop,
												vector w, 
												array[] int v,
												array[] int u,
												int n_cols,
                        vector oset,
                        array[] vector beta_local) {
    real partial_sum = 0.0;


    int n_samples = num_elements(slice_samples);

    for (s_idx in 1:n_samples) {
			int s = slice_samples[s_idx];
			int y_start = y_start_stop[s,1];
      int y_stop = y_start_stop[s,2];

			int n = y_stop - y_start + 2;
  		array[n] int slice_row_ptr = u[y_start:(y_stop + 1)];
			for (i in 1:n) {
				slice_row_ptr[i] = slice_row_ptr[i] - u[y_start] + 1;
			}

			int data_start = data_start_stop[s,1];
      int data_stop = data_start_stop[s,2];

			int n_rows = n - 1;
			
      // matrix vector multiplication
      vector[n_rows] Xb = csr_matrix_times_vector(n_rows,n_cols,w[data_start:data_stop],v[data_start:data_stop],slice_row_ptr,beta_local[s]);

      partial_sum += bernoulli_logit_lupmf(is_cell[y_start:y_stop] | Xb + oset[y_start:y_stop]);
    }

    return partial_sum;
  }
}


data {
  // Basics
  int<lower=1> num_indiv;              // number of patients
  int<lower=1> num_types;               // number of cell types (or interaction targets)
  int<lower=1> num_pot;                 // number of potentials
  int<lower=0> num_pt_groups;           // number of patient groups (can be 0)
  int<lower=0> n_cells;                 // number of cells
  int<lower=0> d_cells;                 // number of spatial features (input dimension)

  // Outcome
  array[n_cells] int<lower=0,upper=1> is_cell;
  vector[n_cells] oset;

  // Sample structure
  int<lower=1> n_samples;
  array[n_samples] int<lower=1,upper=num_indiv> sample_to_indiv;
  array[n_samples,2] int<lower=1,upper=n_cells> y_start_stop;

  // Sparse matrix (CSR) structure
  int n_nz;
  array[n_samples,2] int<lower=1,upper=n_nz> data_start_stop;
  vector[n_nz] w;
  array[n_nz] int<lower=1,upper=n_cells> v;
  array[n_cells + 1] int<lower=1,upper=n_nz+1> u;

  // Patient structure
  array[num_indiv] int<lower=0,upper=num_pt_groups> indiv_to_group;

  // Hyperparameters
  real mean_alpha;
  vector<lower=0>[num_pot] scale_sigma_betas;
  real<lower=0> scale_sigma_alpha;
  real<lower=0> scale_sigmas;

  // Parallelization
  int<lower=1> grainsize;
}

transformed data {

  int num_combos = num_types - 1;
  array[num_combos,num_pot] int beta_idx;
  int beta_start = 1;
  for(i in 1:num_combos) {
    for(j in 1:num_pot) {
      beta_idx[i,j] = (i-1) * num_pot + j + beta_start;
    }
  }

  array[n_samples] int samples_list;
  for(s in 1:n_samples) {
    samples_list[s] = s;
  }
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

}

model {
  // --- Patient priors ---
  if (num_pt_groups > 0) {
    for (j in 1:num_pot) {
      sigma_beta_global[j] ~ normal(scale_sigma_betas[j], scale_sigma_betas[j]);
      to_vector(beta_global[beta_idx[:,j], :]) ~ normal(0, sigma_beta_global[j]);
    }
    tau_alpha_global ~ normal(scale_sigma_alpha, 10);
    to_vector(beta_global[1, :]) ~ normal(mean_alpha, tau_alpha_global[1]);

    sigma_beta_indiv ~ normal(0, scale_sigmas);
    for (i in 1:num_indiv) {
      beta_indiv[:,i] ~ normal(beta_global[, indiv_to_group[i]], sigma_beta_indiv[1]);
    }
  } else {
    for (j in 1:num_pot) {
      sigma_beta_indiv[j] ~ normal(scale_sigma_betas[j], scale_sigma_betas[j]);
      to_vector(beta_indiv[beta_idx[:,j], :]) ~ normal(0, sigma_beta_indiv[j]);
    }
    to_vector(beta_indiv[1, :]) ~ normal(mean_alpha, tau_alpha_indiv);
  }

  // --- Sample priors ---
  sigma_beta_local ~ normal(0, scale_sigmas);
  tau_alpha_indiv ~ normal(0, scale_sigmas);

  for (s in 1:n_samples) {
    beta_local[s] ~ normal(
      beta_indiv[:, sample_to_indiv[s]],
      sigma_beta_local
    );
  }

  // --- Likelihood ---
  target += reduce_sum(partial_sum_lpmf, samples_list, grainsize,
                       is_cell, y_start_stop, data_start_stop,
                       w, v, u, d_cells,
                       oset, beta_local);
}
