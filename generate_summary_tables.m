% GENERATE_SUMMARY_TABLES  Post-processing of the results_*.mat files.
%
% Run once after the case scripts (case33, case34_conical_island for A/B/C,
% case35_inclined_bay) have written their results_<case>.mat files. Prints
%   Table 1: error metrics against reference / analytical data
%   Table 2: mass conservation and positivity
%   Table 3: computational cost
%
% Requires: compute_error_metrics.m on the MATLAB path;
%           Data_NOAAori.xlsx, dataB.txt, dataC.txt, cg_eksak.mat

case_files = { ...
    'results_carrier_greenspan.mat', ...
    'results_conical_island_A.mat', ...
    'results_conical_island_B.mat', ...
    'results_conical_island_C.mat', ...
    'results_inclined_bay.mat' };

%% Reference data
% Conical island: gauge time series scaled to the dimensionless time
% t*sqrt(g/depth) and elevation eta/depth.
depth   = 0.32;                % still-water depth of the experiment [m]
g_const = 9.81;

% Case A (NOAA data, sheet TS2A)
G_A = readtable('Data_NOAAori.xlsx','Sheet','TS2A','Range','A6:I1506');
t_ref_A   = G_A.Time*sqrt(g_const/depth);
eta_ref_A = [G_A.g1_m, G_A.g2_m, G_A.g3_m, G_A.g4_m, ...
             G_A.g6_m, G_A.g9_m, G_A.g16_m, G_A.g22_m]/depth;
ref.conical_island_A = arrayfun(@(gi) {t_ref_A, eta_ref_A(:,gi)}, 1:8, 'UniformOutput', false);

% Case B and C: 9 columns without header = time, WG1,2,3,4,6,9,16,22
D_B = readmatrix('dataB.txt');
ref.conical_island_B = arrayfun(@(gi) {D_B(:,1)*sqrt(g_const/depth), D_B(:,gi+1)/depth}, ...
                                1:8, 'UniformOutput', false);
D_C = readmatrix('dataC.txt');
ref.conical_island_C = arrayfun(@(gi) {D_C(:,1)*sqrt(g_const/depth), D_C(:,gi+1)/depth}, ...
                                1:8, 'UniformOutput', false);

% Carrier-Greenspan: analytical solution, columns = (t, eta)
S_cg      = load('cg_eksak.mat');
data_cg   = S_cg.A;
ref.carrier_greenspan = {{data_cg(:,1), data_cg(:,2)}};

% Inclined bay: no reference data loaded yet. Reference = digitized
% run-up curves of Hubbard & Dodd (2002), Fig. 9, as a 1x5 cell of
% {t_ref, eta_ref} pairs in raw (not normalized) units.
ref.inclined_bay = [];

%% Time shift between simulation clock and reference clock
delt.conical_island_A = 148.59 - 4.9;
delt.conical_island_B = 144.22;
delt.conical_island_C = 144.5;
delt.carrier_greenspan = 0;
delt.inclined_bay      = 0;

%% Comparison window [t_lo t_hi] (dimensionless time) for the error metrics
comp_window.conical_island_A = [144 221.7];
comp_window.conical_island_B = [144 221.7];
comp_window.conical_island_C = [144 221.7];
comp_window.carrier_greenspan = [-Inf Inf];
comp_window.inclined_bay      = [-Inf Inf];

%% Table 1: error metrics
fprintf('=== TABLE 1: Error metrics ===\n');
fprintf('%-22s %-8s %10s %10s %8s\n', 'Case', 'Gauge#', 'RMSE', 'NRMSE[%]', 'R');
for k = 1:numel(case_files)
    if ~isfile(case_files{k})
        fprintf('  [skip] %s not found yet\n', case_files{k});
        continue;
    end
    R = load(case_files{k});
    case_ref = ref.(R.case_name);
    if isempty(case_ref)
        fprintf('  [skip] no reference data loaded yet for %s -- fill in ref.%s above\n', ...
                R.case_name, R.case_name);
        continue;
    end
    for g = 1:size(R.eta_gauge, 2)
        t_ref_full = case_ref{g}{1};
        y_ref_full = case_ref{g}{2};

        % crop the reference series to the comparison window
        w = comp_window.(R.case_name);
        in_window = t_ref_full >= w(1) & t_ref_full <= w(2);
        t_ref = t_ref_full(in_window);
        y_ref = y_ref_full(in_window);

        % Simulation time -> dimensionless time of the reference, then shift.
        % (cg_eksak.mat is already on the same time scale as R.t_hist.)
        if strcmp(R.case_name, 'carrier_greenspan')
            t_sim = R.t_hist + delt.(R.case_name);
        else
            t_sim = R.t_hist*sqrt(g_const) + delt.(R.case_name);
        end

        % warn if the simulation does not reach the end of the window
        target_hi = min(w(2), max(t_ref_full));
        if isfinite(target_hi) && max(t_sim) < target_hi
            fprintf('  [WARNING] %s gauge %d: t_sim only reaches %.2f, window ends at %.2f -- error computed on a shorter window\n', ...
                    R.case_name, g, max(t_sim), target_hi);
        end

        m = compute_error_metrics(t_sim, R.eta_gauge(:,g), t_ref, y_ref);
        fprintf('%-22s %-8d %10.4e %10.2f %8.3f\n', ...
                R.case_name, g, m.RMSE, m.NRMSE, m.R);
    end
end

%% Table 2: conservation and positivity
% The history arrays have Nt+1 entries but the time loop fills only
% 1..Nt, so the last entry (still 0) is dropped.
fprintf('\n=== TABLE 2: Conservation and positivity ===\n');
fprintf('%-22s %14s %14s %16s\n', 'Case', 'MaxVolDrift[%]', 'MinH_after[m]', 'ClipRate[%]');
for k = 1:numel(case_files)
    if ~isfile(case_files{k}); continue; end
    R = load(case_files{k});
    Vh = R.Volume_hist(1:R.Nt);
    volDrift = max(abs((Vh - Vh(1))/Vh(1)*100));
    clipRate = 100*sum(R.nClipped_hist(1:R.Nt))/(R.Nx*R.Ny*R.Nt);
    fprintf('%-22s %14.3e %14.3e %16.4f\n', ...
            R.case_name, volDrift, min(R.minH_hist(1:R.Nt)), clipRate);
end

%% Table 3: computational cost
% Cases A/B/C share the same grid and dt, but their per-step cost can
% differ, so each case is reported with its own numbers.
fprintf('\n=== TABLE 3: Computational cost ===\n');
fprintf('%-22s %10s %8s %12s %12s %10s\n', ...
        'Case', 'Grid', 'Nt', 'Total[s]', 's/step', 'Fact.share');
for k = 1:numel(case_files)
    if ~isfile(case_files{k}); continue; end
    R = load(case_files{k});
    factShare = mean(R.t_fact_hist)/(mean(R.t_fact_hist)+mean(R.t_solve_hist));
    fprintf('%-22s %4dx%-5d %8d %12.2f %12.4f %10.1f%%\n', ...
            R.case_name, R.Nx, R.Ny, R.Nt, R.T_total, R.T_total/R.Nt, 100*factShare);
end
