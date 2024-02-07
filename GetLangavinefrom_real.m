%% Explore the dynamics in terms of alternative stable states and tipping point


for i = 1:806
try
% Get the sleep duration values 
load('Complied_sleep_timesJan2024.mat', ['TapSub_', num2str(i)])
x = eval(['TapSub_', num2str(i)]);


% Get daily values 
x_ = diff(single(getdailysleepdur(x, 8)));
if length(x_) > 100
% Fit langavine dynamics 


    % z-scoring the trajectories - allows us to come back to the age scale
    % later
    MU = mean(x_, 'omitnan');
    SIGMA = std(x_, 'omitnan');
    x_ = (x_ - MU) / SIGMA;
    all_MU = MU;
    all_SIGMA = SIGMA;

    avec = linspace(min(x_), max(x_), length(x_) / 5);
    bw = 0.3 * std(x_, 'omitnan');
    DT = 1;
    nx = length(diff(x_));
    % Calculate the Langevin reconstruction by means of the mesh method 
    % as opposed to the bin method in Arani et al.
    results_M = LangevinReconst_MESH(x_,diff(x_),nx,DT,bw,length(avec),avec);

    % get the model from the results
    mod = langevin_eq(results_M);
mod.plot('potential_eff');
uiwait
    % use functionality to calculate potential eff
    %ueff = mod.potential_eff.ueff;
    %U = ueff(mod.potential_eff.dom);
    %normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));
end
end
end
    % 
    % % sometimes we get complex numbers, model has failed, we skip sample
    % if ~isreal(U)
    %     continue
    % end
    % 
    % % first order approximation of the average derivative of the DT landscape
    % derivative_x =  diff(U);
    % derivative_patterns(i) = mean(derivative_x);
    % 
    % all_dom_ranges(i) = mod.potential_eff.dom(end) - mod.potential_eff.dom(1);
    % 
    % % extract depth - width - slope
    % % calculate equilibria of the effective potential landascape
    % mod.equilibria = mod.find_equilibria('effective');
    % 
    % n_equilibria(i) = length(mod.equilibria);
    % 
    % % extract stable and unstable points
    % unstable_eq = mod.equilibria([mod.equilibria.stable] == 0);
    % stable_eq = mod.equilibria([mod.equilibria.stable] == 1);
    % 
    % all_stable_points{i} = stable_eq;
    % 
    % n_stable_points(i) = length(stable_eq);
    % n_unstable_points(i) = length(unstable_eq);
    % 
    % % wheter or not an unstable point is present
    % all_has_unstable(i) = ~isempty(unstable_eq);
    % 
    % % pre-allocation
    % slopes = zeros(1, length(stable_eq));
    % eq = cell(1, length(stable_eq));
