% https://github.com/libDirectional/libDirectional 

clear;
rng(8675309);

% =========================================================================
% load data
% =========================================================================

data_degrees = readmatrix('datasets/black_mountain_wind_direction.csv');
data_radians = data_degrees * pi / 180;
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
Re2 = (T0 / (T0 - 1)) * (Rbar^2 - 1/T0);
sig = sqrt(log(1/Re2));
kappa = Rbar * (n - Rbar^2) / (1 - Rbar^2);

% =========================================================================
% pre-allocate storage
% =========================================================================

vmf_forecast_draws = zeros(T, M);
wn_forecast_draws = zeros(T, M);

% =========================================================================
% initialize filters
% =========================================================================

vmf_filter = VMFFilter();
vmf_filter.setState(VMFDistribution([1; 0], 1));

wn_filter = WNFilter();
wn_filter.setState(WNDistribution(0, 1));

for t = 1:T

  % =======================================================================
  % compute prior predictive distribution of state
  % =======================================================================

  vmf_filter.predictIdentity(VMFDistribution([0; 1], 1));
  wn_filter.predictIdentity(WNDistribution(0, 1))

  vmf_prior = vmf_filter.getEstimate();
  wn_prior = wn_filter.getEstimate();

  % =======================================================================
  % simulate one-step-ahead predictive distribution
  % =======================================================================

  for m = 1:M

    vmf_s_draw = sample(vmf_prior, 1);
    vmf_y_draw = sample(VMFDistribution(vmf_s_draw, kappa), 1);

    wn_s_draw = sample(wn_prior, 1);
    wn_y_draw = mod(wn_s_draw + sample(WNDistribution(0, sig), 1), 2 * pi);

    vmf_forecast_draws(t, m) = mod(atan2(vmf_y_draw(2), vmf_y_draw(1)), 2 * pi);
    wn_forecast_draws(t, m) = wn_y_draw;

  end

  % =======================================================================
  % compute filtering distribution
  % =======================================================================

  vmf_filter.updateIdentity(VMFDistribution([0; 1], kappa), data_unitvec(t, :)');
  wn_filter.updateIdentity(WNDistribution(0, sig), data_radians(t))

end

% =========================================================================
% preallocate storage
% =========================================================================

writematrix(vmf_forecast_draws, 'from_matlab_vmf_black_mountain_forecasts.csv')
writematrix(wn_forecast_draws, 'from_matlab_wn_black_mountain_forecasts.csv')
