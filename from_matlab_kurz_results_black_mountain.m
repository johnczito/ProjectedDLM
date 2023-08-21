% https://github.com/libDirectional/libDirectional 

data_radians = readmatrix('black_mountain_radians.csv');
data_unitvec = [cos(data_radians) sin(data_radians)];
T = length(data_radians);
M = 1000;
vmf_forecast_draws = zeros(T, M);
wn_forecast_draws = zeros(T, M);

vmf_filter = VMFFilter();
vmf_filter.setState(VMFDistribution([1; 0], 1));

wn_filter = WNFilter();
wn_filter.setState(WNDistribution(0, 1));

for t = 1:T

  vmf_filter.predictIdentity(VMFDistribution([0; 1], 1));
  wn_filter.predictIdentity(WNDistribution(0, 1))

  vmf_prior = vmf_filter.getEstimate();
  wn_prior = wn_filter.getEstimate();

  for m = 1:M

    vmf_s_draw = sample(vmf_prior, 1);
    vmf_y_draw = sample(VMFDistribution(vmf_s_draw, 1), 1);

    wn_s_draw = sample(wn_prior, 1);
    wn_y_draw = mod(wn_s_draw + sample(WNDistribution(0, 1), 1), 2 * pi);

    vmf_forecast_draws(t, m) = mod(atan2(vmf_y_draw(2), vmf_y_draw(1)), 2 * pi);
    wn_forecast_draws(t, m) = wn_y_draw;

  end

  vmf_filter.updateIdentity(VMFDistribution([0; 1], 1), data_unitvec(t, :)');
  wn_filter.updateIdentity(WNDistribution(0, 1), data_radians(t))

end

writematrix(vmf_forecast_draws, 'from_matlab_vmf_black_mountain_forecasts.csv')
writematrix(wn_forecast_draws, 'from_matlab_wn_black_mountain_forecasts.csv')
