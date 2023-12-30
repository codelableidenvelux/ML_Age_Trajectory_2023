%% null distribution
% Following https://www.pnas.org/doi/10.1073/pnas.0802430105#supplementary-materials
% We want to see in what cases a surrgate time series presents alternative
% stable states by chance. 
% for each sample:
    % 1. We create 1000 surrogate times series (each of 3 methods)
    % 2. We calculate the potential_eff
    % 3. We see what percentage of those have >= number of stable points
    % w.r.t the original time series
    % 4. We see in how many of the samples the results has P. 0.05
% Now for the 3 methods
% a. simple randomization with replacement 
% b. same spectrum (?)
% c. autoregressive model

load('all_age_pred_21d_30_11_2023.mat')

[sorted_ages, idxs] = sort(cellfun(@(x) x(end), all_ages));
sorted_genders = all_genders(idxs);
sorted_preds = all_preds(idxs);
sorted_pred_times = all_preds_times(idxs);

sorted_pids = all_pids(idxs);

%% A. randomization
n_sub = length(all_preds);
n_rep = 1000;
surrogate_models_rnd = cell(n_rep, n_sub);

for i = 1:n_sub

    % get trajectory for subject
    x_ = double(sorted_preds{i});
    % z-scoring the trajectories 
    MU = mean(x_, 'omitnan');
    SIGMA = std(x_, 'omitnan');
    x_ = (x_ - MU) / SIGMA;

    ys_ = surrogate_rnd(x_, n_rep);

    for j = 1:n_rep

        avec = linspace(min(ys_(j, :)), max(ys_(j, :)), length(ys_(j, :)) / 5);
        bw = 0.3 * std(ys_(j, :), 'omitnan');
        DT = 1;
        nx = length(diff(ys_(j, :)));

        % Calculate the Langevin reconstruction by means of the mesh method 
        % as opposed to the bin method in Arani et al.
        results_M = LangevinReconst_MESH(ys_(j, :),diff(ys_(j, :)),nx,DT,bw,length(avec),avec);
    
        % get the model from the results
        mod = langevin_eq(results_M);
        surrogate_models_rnd{j, i} = mod;
    end

end


%% C. Autoregressive model
addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
n_sub = length(all_preds);
n_rep = 1000;

for i = 241:n_sub
    surrogate_models_ar = cell(n_rep, 1);
    % get trajectory for subject
    x_ = double(sorted_preds{i});
    % z-scoring the trajectories 
    MU = mean(x_, 'omitnan');
    SIGMA = std(x_, 'omitnan');
    x_ = (x_ - MU) / SIGMA;

    avec = linspace(min(x_), max(x_), length(x_) / 5);
    bw = 0.3 * std(x_, 'omitnan');
    DT = 1;
    nx = length(diff(x_));
    % Calculate the Langevin reconstruction by means of the mesh method 
    % as opposed to the bin method in Arani et al.
    results_M = LangevinReconst_MESH(x_,diff(x_),nx,DT,bw,length(avec),avec);

    % get the model from the results
    mod = langevin_eq(results_M);

    ys_ = surrogate_autoreg(x_, n_rep);

    parfor j = 1:n_rep
        addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
        fprintf("Rep %d/%d || SUB[%d/%d]\n", j, n_rep, i, n_sub)
        avec = linspace(min(ys_(j, :)), max(ys_(j, :)), length(ys_(j, :)) / 5);
        bw = 0.3 * std(ys_(j, :), 'omitnan');
        DT = 1;
        nx = length(diff(ys_(j, :)));

        % Calculate the Langevin reconstruction by means of the mesh method 
        % as opposed to the bin method in Arani et al.
        results_M = LangevinReconst_MESH(ys_(j, :),diff(ys_(j, :)),nx,DT,bw,length(avec),avec);
        
        % get the model from the results
        mod = langevin_eq(results_M);
        mod.equilibria = mod.find_equilibria('effective');

        if (~isempty(mod.equilibria([mod.equilibria.stable] == 0)))
            mean_exit_all = mod.mean_exit('all', mod.pdf);
        end

        surrogate_models_ar{j} = mod;
    end

    save(sprintf('./surrogate_ar/surrogate_ar_sorted_%d.mat', i), 'mod', 'x_', 'ys_', 'MU', 'SIGMA', 'surrogate_models_ar')

end 

%% Left most stable point
all_pctile = [];
for idx = 1:144
    try
        load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ar/surrogate_ar_sorted_', num2str(idx), '.mat'])
        mod = all_mods{idx};
        pos_left = zeros(1000, 1);
        
        for i = 1:1000
            stable_points = surrogate_models_ar{i}.equilibria([surrogate_models_ar{i}.equilibria.stable] == 1);
            pos_left(i) = stable_points(1).x;
        end
        
        stable_points = mod.equilibria([mod.equilibria.stable] == 1);
        
        
        % histogram(pos_left, 'Normalization','probability')
        % hold on
        % plot([ stable_points(1).x, stable_points(1).x ], [0, .5], '--k')
        % plot([ prctile(pos_left, 95), prctile(pos_left, 95) ], [0, .5], '--r')
        % 
        % Test data
        historicalData = rand(1000, 1);
        exogenousVariable = 0.7;
        
        % Compute centile
        nless = sum(pos_left < stable_points(1).x);
        nequal = sum(pos_left == stable_points(1).x);
        centile = 100 * (nless + 0.5*nequal) / length(pos_left);
        all_pctile = [all_pctile, centile];
    end
end

%% Left most UNstable point
all_pctile_unstable = [];
for idx = 1:144
    try
        load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ar/surrogate_ar_sorted_', num2str(idx), '.mat'])
        fprintf('%d\n', idx)
        mod = all_mods{idx};
        unstable_points = mod.equilibria([mod.equilibria.stable] == 0);
        if length(unstable_points) >= 1
            pos_left = [];

            for i = 1:1000
                unstable_points = surrogate_models_ar{i}.equilibria([surrogate_models_ar{i}.equilibria.stable] == 0);
                if length(unstable_points) >= 1
                    pos_left = [pos_left, unstable_points(1).x];
                else
                    pos_left = [pos_left, -Inf];
                end
            end

            unstable_points = mod.equilibria([mod.equilibria.stable] == 0);

            % Compute centile
            nless = sum(pos_left < unstable_points(1).x);
            nequal = sum(pos_left == unstable_points(1).x);
            centile = 100 * (nless + 0.5*nequal) / length(pos_left);
            all_pctile_unstable = [all_pctile_unstable, centile];
        end
    catch ME
        
    end
end
%% max(stable points difference)
all_pctile_diff = nan(164, 1);
for idx = 1:164
    try
        load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ar/surrogate_ar_sorted_', num2str(idx), '.mat'])
        fprintf('%d\n', idx)
        stable_points = mod.equilibria([mod.equilibria.stable] == 1);
        % we only care about people with at least 2 stable points
        if length(stable_points) > 1
            pos_left = [];

            for i = 1:1000
                stable_points_sur = surrogate_models_ar{i}.equilibria([surrogate_models_ar{i}.equilibria.stable] == 1);
                if length(stable_points) == length(stable_points_sur)
                    pos_left = [ pos_left, abs(stable_points(end).x - stable_points(1).x)];
                else
                    pos_left = [ pos_left, 0];
                end

            end
            stable_points = mod.equilibria([mod.equilibria.stable] == 1);
            real_diff = abs(stable_points(end).x - stable_points(1).x);

            % Compute centile
            nless = sum(pos_left < real_diff);
            nequal = sum(pos_left == real_diff);
            centile = 100 * (nless + 0.5*nequal) / length(pos_left);
            all_pctile_diff(idx) = centile;
            fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)

        end
    end
end

%% max(MET)
all_pctile_met = nan(164, 1);
for idx = 1:164
    try
        load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ar/surrogate_ar_sorted_', num2str(idx), '.mat'])
        mod = all_mods{idx};
        stable_points = mod.equilibria([mod.equilibria.stable] == 1);
        if length(stable_points) > 1
            pos_left = [];

            for i = 1:1000
                stable_points_sur = surrogate_models_ar{i}.equilibria([surrogate_models_ar{i}.equilibria.stable] == 1);
                
                if length(stable_points) == length(stable_points_sur)

                    mean_exit_all = surrogate_models_ar{i}.xtra.results.mean_exit;
                    val = max(cellfun(@(x) x.WT, mean_exit_all));

                    pos_left = [pos_left, val];
                else
                    pos_left = [pos_left, 0];
                end
            end

            
            mean_exit_all = mod.mean_exit('all', mod.pdf);
            val = max(cellfun(@(x) x.WT, mean_exit_all));
            % Compute centile
            nless = sum(pos_left < val);
            nequal = sum(pos_left == val);
            centile = 100 * (nless + 0.5*nequal) / length(pos_left);
            all_pctile_met(idx) = centile;
            fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)
        end
    catch ME
        
    end
end