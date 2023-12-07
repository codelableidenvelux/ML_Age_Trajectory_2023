x0 = sorted_preds{123};

MU = mean(x0, 'omitnan');
SIGMA = std(x0, 'omitnan');
x1 = double((x0 - MU) / SIGMA);

% n_bins = 30;
% 
% results = LangevinReconst(x1, prctile(x1,1), prctile(x1,99), n_bins, 1:5, 1, 'Nadaraya-Watson');

avec = linspace(min(x1), max(x1), length(x1) / 5);
bw = 0.3 * std(x1);
DT = 1;
nx = length(diff(x1));
results_M = LangevinReconst_MESH(x1,diff(x1),nx,DT,bw,length(avec),avec);
mod_M = langevin_eq(results_M);
subplot(2,1,1)
plot(x0)
subplot(2,1,2)
mod_M.plot('potential_eff')
%%
for i = 1:100

    x2 = x1(randperm(length(x1)));
    results_S{i} = LangevinReconst_MESH(x2,diff(x2),nx,DT,bw,length(avec),avec);
    
end

results_S = LangevinReconst_MESH(x1,diff(x1),nx,DT,bw,length(avec),avec);

mod_M = langevin_eq(results_M);
mod = langevin_eq(results);

figure(5)
mod.plot('potential_eff')