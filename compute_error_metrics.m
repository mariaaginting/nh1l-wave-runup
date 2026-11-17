function metrics = compute_error_metrics(t_sim, y_sim, t_ref, y_ref)
%COMPUTE_ERROR_METRICS  RMSE, NRMSE and correlation between a simulated and
%a reference (experimental/analytical) time series.
%
%   metrics = compute_error_metrics(t_sim, y_sim, t_ref, y_ref)
%
%   INPUTS
%     t_sim, y_sim : time vector and simulated signal (e.g. eta(t) at one gauge)
%     t_ref, y_ref : time vector and reference signal. t_ref need not share
%                    the sampling of t_sim: y_sim is interpolated onto t_ref
%                    (reference points outside the simulated time span are
%                    dropped).
%
%   OUTPUT  (struct)
%     metrics.RMSE    root-mean-square error, same units as y
%     metrics.NRMSE   RMSE / (max(y_ref)-min(y_ref)) * 100   [%]
%     metrics.R       Pearson correlation coefficient
%     metrics.N       number of points used
%     metrics.t, metrics.y_sim_i, metrics.y_ref   (for plotting)

t_ref = t_ref(:);
y_ref = y_ref(:);
t_sim = t_sim(:);
y_sim = y_sim(:);

valid = t_ref >= min(t_sim) & t_ref <= max(t_sim);
t_ref = t_ref(valid);
y_ref = y_ref(valid);

y_sim_i = interp1(t_sim, y_sim, t_ref, 'linear');

err = y_sim_i - y_ref;
N   = numel(err);

RMSE = sqrt(mean(err.^2));
rng  = max(y_ref) - min(y_ref);
if rng > 0
    NRMSE = RMSE/rng*100;
else
    NRMSE = NaN;
end

ys = y_sim_i - mean(y_sim_i);
yr = y_ref   - mean(y_ref);
R  = sum(ys.*yr)/sqrt(sum(ys.^2)*sum(yr.^2));

metrics.RMSE    = RMSE;
metrics.NRMSE   = NRMSE;
metrics.R       = R;
metrics.N       = N;
metrics.t       = t_ref;
metrics.y_sim_i = y_sim_i;
metrics.y_ref   = y_ref;
end
