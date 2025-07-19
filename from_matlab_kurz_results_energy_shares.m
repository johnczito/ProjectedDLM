% https://github.com/libDirectional/libDirectional 

clear;
rng(8675309);

% =========================================================================
% load data
% =========================================================================

which_dim = "high"

if which_dim == "low"
  data_unitvec = readmatrix('datasets/annual_3d_shares_1990_2022.csv');
else 
  data_unitvec = readmatrix('datasets/annual_8d_shares_1990_2022.csv');
end

% =========================================================================
% store dimensions
% =========================================================================

T = size(data_unitvec, 1);
n = size(data_unitvec, 2);
M = 2500;

% =========================================================================
% estimate static parameters on a pre-sample
% =========================================================================

T0 = 10;

z = [zeros(n-1,1); 1];

Rbar = norm(mean(data_unitvec(1:T0, :), 1));

kappa = Rbar * (n - Rbar^2) / (1 - Rbar^2);

% =========================================================================
% pre-allocate storage
% =========================================================================

vmf_forecast_draws = zeros(M, n, T);

% =========================================================================
% initialize filters
% =========================================================================

vmf_filter = VMFFilter();
vmf_filter.setState(VMFDistribution(z, 1));

for t = 1:T

  vmf_filter.predictIdentity(VMFDistribution(z, 1));

  vmf_prior = vmf_filter.getEstimate();

  for m = 1:M

    vmf_s_draw = sample(vmf_prior, 1);
    vmf_y_draw = sample(VMFDistribution(vmf_s_draw, kappa), 1);

    vmf_forecast_draws(m, :, t) = vmf_y_draw;

  end

  vmf_filter.updateIdentity(VMFDistribution(z, kappa), data_unitvec(t, :)');

end

save("from_matlab_kurz_vmf_energy_shares_forecasts.mat", 'vmf_forecast_draws')
