% https://github.com/libDirectional/libDirectional 

clear;
rng(8675309);

% =========================================================================
% load data
% =========================================================================

data_radians = readmatrix('datasets/ohare_wind_direction.csv');
data_unitvec = [cos(data_radians) sin(data_radians)];

% =========================================================================
% store dimensions
% =========================================================================

T = length(data_radians);
n = size(data_unitvec, 2);
M = 2500;

% =========================================================================
% estimate static parameters on a pre-sample
% =========================================================================

T0 = T;

Rbar = norm(mean(data_unitvec(1:T0, :), 1));

kappa = Rbar * (n - Rbar^2) / (1 - Rbar^2);

% =========================================================================
% pre-allocate storage
% =========================================================================

vmf_forecast_draws = zeros(T, M);

% =========================================================================
% initialize filters
% =========================================================================

vmf_filter = VMFFilter();
vmf_filter.setState(VMFDistribution([1; 0], 1));

for t = 1:T

  % =======================================================================
  % compute prior predictive distribution of state
  % =======================================================================

  vmf_filter.predictIdentity(VMFDistribution([0; 1], 1));

  vmf_prior = vmf_filter.getEstimate();

  % =======================================================================
  % simulate one-step-ahead predictive distribution
  % =======================================================================

  for m = 1:M

    vmf_s_draw = sample(vmf_prior, 1);
    vmf_y_draw = sample(VMFDistribution(vmf_s_draw, kappa), 1);

    vmf_forecast_draws(t, m) = mod(atan2(vmf_y_draw(2), vmf_y_draw(1)), 2 * pi);

  end

  % =======================================================================
  % compute filtering distribution
  % =======================================================================

  vmf_filter.updateIdentity(VMFDistribution([0; 1], kappa), data_unitvec(t, :)');

end

% =========================================================================
% preallocate storage
% =========================================================================

writematrix(vmf_forecast_draws, 'from_matlab_vmf_ohare_forecasts.csv')
