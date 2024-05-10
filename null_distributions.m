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
% b. same spectrum
% c. autoregressive model

load('all_age_pred_21d_30_11_2023.mat')
load('analysis_output_landscapes_15_12_2023.mat', 'all_mods')

[sorted_ages, idxs] = sort(cellfun(@(x) x(end), all_ages));
sorted_genders = all_genders(idxs);
sorted_preds = all_preds(idxs);
sorted_pred_times = all_preds_times(idxs);

sorted_pids = all_pids(idxs);

%% A. randomization
addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
n_sub = length(all_preds);
n_rep = 1000;

for i = 1:n_sub
    surrogate_models = cell(n_rep, 1);
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
    mod.equilibria = mod.find_equilibria('effective');
    stable_points = mod.equilibria([mod.equilibria.stable] == 1);

    ys_ = surrogate_rnd(x_, n_rep);

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
        mod_ = langevin_eq(results_M);
        mod_.equilibria = mod_.find_equilibria('effective');

        stable_points_sur = mod_.equilibria([mod_.equilibria.stable] == 1);
        try
            if length(stable_points) == length(stable_points_sur)
                if (~isempty(mod_.equilibria([mod_.equilibria.stable] == 0)))
                    mean_exit_all = mod_.mean_exit('all', mod_.pdf);
                end
            end
        catch err
            disp(err.message);
        end

        surrogate_models{j} = mod_;
    end
        save(sprintf('./surrogate_rnd/surrogate_rnd_sorted_%d.mat', i), 'mod', 'x_', 'ys_', 'MU', 'SIGMA', 'surrogate_models')
end


%% B. Preserve spectrum (FT)
addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
addpath '/home/enea/code/ML_Age_Trajectory_2023'/LancsBioMed/Surrogate_Data_Testing/SurrogateData/
n_sub = length(all_preds);
n_rep = 1000;

for i = 1:n_sub

    surrogate_models = cell(n_rep, 1);
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
    mod.equilibria = mod.find_equilibria('effective');

    stable_points = mod.equilibria([mod.equilibria.stable] == 1);

    ys_ = surrogate(x_, n_rep, 'FT', 0, 1);

    parfor j = 1:n_rep
        addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
        addpath '/home/enea/code/ML_Age_Trajectory_2023'/LancsBioMed/Surrogate_Data_Testing/SurrogateData/
        fprintf("Rep %d/%d || SUB[%d/%d]\n", j, n_rep, i, n_sub)
        avec = linspace(min(ys_(j, :)), max(ys_(j, :)), length(ys_(j, :)) / 5);
        bw = 0.3 * std(ys_(j, :), 'omitnan');
        DT = 1;
        nx = length(diff(ys_(j, :)));

        % Calculate the Langevin reconstruction by means of the mesh method 
        % as opposed to the bin method in Arani et al.
        results_M = LangevinReconst_MESH(ys_(j, :),diff(ys_(j, :)),nx,DT,bw,length(avec),avec);
    
        % get the model from the results
        mod_ = langevin_eq(results_M);
        mod_.equilibria = mod_.find_equilibria('effective');
        try
            stable_points_sur = mod_.equilibria([mod_.equilibria.stable] == 1);
        
            if length(stable_points) == length(stable_points_sur)
                if (~isempty(mod_.equilibria([mod_.equilibria.stable] == 0)))
                    mean_exit_all = mod_.mean_exit('all', mod_.pdf);
                end
            end
        catch err
            disp(err.message);
        end

        surrogate_models{j} = mod_;
    end
        save(sprintf('./surrogate_ft/surrogate_ft_sorted_%d.mat', i), 'mod', 'x_', 'ys_', 'MU', 'SIGMA', 'surrogate_models')
end


%% C. Autoregressive model
addpath '/home/enea/code/ML_Age_Trajectory_2023/aay4895-matlab-code'
n_sub = length(all_preds);
n_rep = 1000;

for i = 1:n_sub
    try
        surrogate_models = cell(n_rep, 1);
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
            mod_ = langevin_eq(results_M);
            mod_.equilibria = mod_.find_equilibria('effective');

            if (~isempty(mod_.equilibria([mod_.equilibria.stable] == 0)))
                mean_exit_all = mod_.mean_exit('all', mod_.pdf);
            end

            surrogate_models{j} = mod_;
        end

        save(sprintf('./surrogate_ar/surrogate_ar_sorted_%d.mat', i), 'mod', 'x_', 'ys_', 'MU', 'SIGMA', 'surrogate_models')
    end
end

%% max(Mean Exit Time) - RND
n_sub = 291;
n_stable_points_across_boots = nan(n_sub, 1000);
all_pctile_met_rnd = nan(n_sub, 1);
all_pct_same_landscape_met_rnd = nan(n_sub, 1);

for idx = 1:n_sub
    try
    load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_rnd/surrogate_rnd_sorted_', num2str(idx), '.mat'])
    mod = all_mods{idx};
    stable_points = mod.equilibria([mod.equilibria.stable] == 1);
    if length(stable_points) > 1
        pos_left = [];

        for i = 1:1000
            stable_points_sur = surrogate_models{i}.equilibria([surrogate_models{i}.equilibria.stable] == 1);
            n_stable_points_across_boots(idx, i) = length(stable_points_sur);
                    
            if length(stable_points) == length(stable_points_sur)
                
                % Some degenerate cases exists where stable/unstable points
                % are duplicated, the number of stable points ends up being
                % the same, but the solution for mean exit time does not 
                % exist .
                if isfield(surrogate_models{i}.xtra, 'results') && isfield(surrogate_models{i}.xtra.results, 'mean_exit')
                    mean_exit_all = surrogate_models{i}.xtra.results.mean_exit;
                    val = max(cellfun(@(x) x.WT, mean_exit_all));
                    
                    pos_left = [pos_left, val];
                else 
                    pos_left = [pos_left, 0];
                end
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
        all_pctile_met_rnd(idx) = centile;
        all_pct_same_landscape_met_rnd(idx) = mean(pos_left > 0) * 100;
        fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)
    end
    catch
    end
end

%% max(Mean Exit Time) - FT
n_sub = 291;
n_stable_points_across_boots = nan(n_sub, 1000);
all_pctile_met_ft = nan(n_sub, 1);
all_pct_same_landscape_met_ft = nan(n_sub, 1);

for idx = 1:n_sub

    load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ft/surrogate_ft_sorted_', num2str(idx), '.mat'])
    mod = all_mods{idx};
    stable_points = mod.equilibria([mod.equilibria.stable] == 1);
    if length(stable_points) > 1
        pos_left = [];

        for i = 1:1000
            stable_points_sur = surrogate_models{i}.equilibria([surrogate_models{i}.equilibria.stable] == 1);
            n_stable_points_across_boots(idx, i) = length(stable_points_sur);
                    
            if length(stable_points) == length(stable_points_sur)
                
                % Some degenerate cases exists where stable/unstable points
                % are duplicated, the number of stable points ends up being
                % the same, but the solution for mean exit time does not 
                % exist .
                if isfield(surrogate_models{i}.xtra, 'results') && isfield(surrogate_models{i}.xtra.results, 'mean_exit')
                    mean_exit_all = surrogate_models{i}.xtra.results.mean_exit;
                    val = max(cellfun(@(x) x.WT, mean_exit_all));
                    
                    pos_left = [pos_left, val];
                else 
                    pos_left = [pos_left, 0];
                end
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
        all_pctile_met_ft(idx) = centile;
        all_pct_same_landscape_met_ft(idx) = mean(pos_left > 0) * 100;
        fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)
    end

end

%% max(Mean Exit Time) - AR
all_pctile_met_ar = nan(n_sub, 1);
all_pct_same_landscape_met_ar = nan(n_sub, 1);

for idx = 1:n_sub
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
            all_pctile_met_ar(idx) = centile;
            all_pct_same_landscape_met_ar(idx) = mean(pos_left > 0) * 100;
            fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)
        end
    catch ME
        
    end
end

%% example
idx = 37;
load(['/home/enea/code/ML_Age_Trajectory_2023/surrogate_ft/surrogate_ft_sorted_', num2str(idx), '.mat'])
mod = all_mods{idx};
stable_points = mod.equilibria([mod.equilibria.stable] == 1);
if length(stable_points) > 1
    pos_left = [];

    for i = 1:1000
        stable_points_sur = surrogate_models{i}.equilibria([surrogate_models{i}.equilibria.stable] == 1);
        
        if length(stable_points) == length(stable_points_sur)

            % get the model from the results
            surrogate_models{i}.equilibria = surrogate_models{i}.find_equilibria('effective');

            if (~isempty(surrogate_models{i}.equilibria([surrogate_models{i}.equilibria.stable] == 0)))
                mean_exit_all = surrogate_models{i}.xtra.results.mean_exit;
                val = max(cellfun(@(x) x.WT, mean_exit_all));
    
                pos_left = [pos_left, val];
            else
                pos_left = [pos_left, 0];
            end

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
    fprintf('%d: %.2f - %.2f\n', idx, centile, mean(pos_left > 0) * 100)
end

histogram(pos_left, 'Normalization','probability')
hold on
plot([ val, val ], [0, .5], '--k')
plot([ prctile(pos_left, 95), prctile(pos_left, 95) ], [0, .5], '--r')
