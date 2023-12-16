% Code to compare the effective potential using the BIN method and the MESH
% method.

% December 2023
% Enea Ceolini (Leiden University & QuantActions AG)
% Arko Ghosh (Leiden University)

% chose trace
x0 = double(sorted_preds{152});

% z-score
MU = mean(x0, 'omitnan');
SIGMA = std(x0, 'omitnan');
x1 = (x0 - MU) / SIGMA;

% parameters MESH
avec = linspace(min(x1), max(x1), length(x1) / 5);
bw = 0.3 * std(x1, 'omitnan');
DT = 1;
nx = length(diff(x1));

% parameters BIN
n_bins = 30;

% BIN
results = LangevinReconst(x1, min(x1), max(x1), n_bins, 1:5, 1, 'Nadaraya-Watson');

% MESH
results_M = LangevinReconst_MESH(x1,diff(x1),nx,DT,bw,length(avec),avec);
mod_M = langevin_eq(results_M);
mod_M.equilibria = mod_M.find_equilibria('effective');
mod = langevin_eq(results);
mod.equilibria = mod.find_equilibria('effective');

% plot
xx1 = mod_M.potential_eff.dom;
xx2 = mod.potential_eff.dom;

xx = [min(xx1(1), xx2(1)), max(xx1(end), xx2(end))];

subplot(3,1,1)
plot(x1, 1:size(x1,2))
xlim(xx)
set(gca, 'YDir','reverse')
subplot(3,1,2)
mod_M.plot('potential_eff')
xlim(xx)
title('MESH method')
subplot(3,1,3)
mod.plot('potential_eff')
xlim(xx)
title('BIN method')
%% null distribution

for i = 1:100
    x2 = x1(randperm(length(x1)));
    results_S{i} = LangevinReconst_MESH(x2,diff(x2),nx,DT,bw,length(avec),avec);
    
end

mod_S = langevin_eq(results_S{42});

figure(5)
mod_S.plot('potential_eff')