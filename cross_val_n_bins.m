%% load healthy age data 20d
n_days = 20;
% load('all_pred_seen_and_unseen_with_crt_20d.mat')
load('all_pred_seen_and_unseen_with_crt_20d_25_04_2023.mat')

[sorted_ages, idxs] = sort(all_ages);
sorted_genders = all_genders(idxs);
sorted_tests_with_time = tests_with_time(idxs);
sorted_preds = all_preds(idxs);
sorted_pred_times = all_preds_times(idxs);
sorted_was_testing = was_testing(idxs);

bins = 7:1:50;
nn_bins = length(bins);

%% analysis

n_subs = length(all_ages);
derivative_patterns = zeros(nn_bins, n_subs);
median_predicted_ages = zeros(nn_bins, n_subs);
gap_real_age_median_predicted_age = zeros(nn_bins, n_subs);
n_equilibria = zeros(nn_bins, n_subs);
n_stable_points = zeros(nn_bins, n_subs);
n_unstable_points = zeros(nn_bins, n_subs);
all_MU = zeros(nn_bins, n_subs);
all_SIGMA = zeros(nn_bins, n_subs);
all_slopes = cell(nn_bins, n_subs);
all_depths = cell(nn_bins, n_subs);
all_widths = cell(nn_bins, n_subs);
all_dists = cell(nn_bins, n_subs);
all_dom_ranges = zeros(nn_bins, n_subs);
all_eq = cell(nn_bins, n_subs);
all_mods = cell(nn_bins, n_subs);
all_has_unstable = zeros(nn_bins, n_subs);

for k = 1:length(bins)
    for i = 1:n_subs
        fprintf("Bins %d/%d (sub %d/%d)\n", k, nn_bins, i, n_subs)
        % get trajectory for subject
        x_ = double(sorted_preds{i});
        x_(x_ == mode(x_)) = NaN;

        % calculate median predicted age over the trajectory
        median_predicted_ages(k, i) = median(x_, 'omitnan');

        % calculate the gap between the median prediction and the real age
        gap_real_age_median_predicted_age(k, i) = median(x_, 'omitnan') - sorted_ages(i);

        % z-scoring the trajectories - allows us to come back to the age scale
        % later
        % CHECK: !!is this correct?!!
        MU = mean(x_, 'omitnan');
        SIGMA = std(x_, 'omitnan');
        x_ = (x_ - MU) / SIGMA;
        all_MU(k, i) = MU;
        all_SIGMA(k, i) = SIGMA;

        % we fix the number of bins depedning on the length of the trajectory
        n_bins = bins(k);


        % calculate the reconstruction model
        % CHECK: !!is this correct?!!
        results = LangevinReconst(x_, prctile(x_,1), prctile(x_,99), n_bins, 1:5, 1, 'Nadaraya-Watson');

        % get the model from the results
        % CHECK: !!is this correct?!!
        mod = langevin_eq(results, 'weightedspline');
        all_mods{k, i} = mod;

        % use functionality to calculate potential eff
        % CHECK: !!is this correct?!!
        ueff = mod.potential_eff.ueff;
        U = ueff(mod.potential_eff.dom);
        normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));
        % sometimes we get complex numbers, does this mean the estimation is
        % wrong?
        % CHECK: !!is this correct?!!
        if ~isreal(U)
            continue
        end

        % first order approximation of the average derivative of the landscape
        derivative_x =  diff(U);
        derivative_patterns(k, i) = mean(derivative_x);

        all_dom_ranges(k, i) = mod.potential_eff.dom(end) - mod.potential_eff.dom(1);

        % extract depth - width - slope
        % calculate equilibria of the effective potential landascape
        mod.equilibria = mod.find_equilibria('effective');

        n_equilibria(k, i) = length(mod.equilibria);

        % extract stable and unstable points
        unstable_eq = mod.equilibria([mod.equilibria.stable] == 0);
        stable_eq = mod.equilibria([mod.equilibria.stable] == 1);

        n_stable_points(k, i) = length(stable_eq);
        n_unstable_points(k, i) = length(unstable_eq);

        % wheter or not an unstable point is present
        all_has_unstable(k, i) = ~isempty(unstable_eq);

        slopes = zeros(1, length(stable_eq));
        depths = zeros(1, length(stable_eq));
        widths = zeros(1, length(stable_eq));
        eq = cell(1, length(stable_eq));

        % here calculate anyway the params even if no unstable points are
        % found.
        if (~isempty(stable_eq)) && (isempty(unstable_eq))

            u_border = [mod.potential_eff.dom(2), mod.potential_eff.dom(end-1)];

            x_eq = stable_eq(1).x;
            y_eq = normU(x_eq);
            [~, idx] = min(abs(u_border - x_eq));
            closest_point = u_border(idx);
            width = abs(x_eq - closest_point);
            y_closest = normU(closest_point);
            depth = y_closest - y_eq;
            slope = depth / width;
            all_widths{k, i} = [width];
            all_depths{k, i} = [depth];
            all_slopes{k, i} = [slope];
            all_eq{k, i} = {[x_eq, y_eq]};
            continue

        end

        for j = 1:length(stable_eq)
            % find closest eq point
            x_eq = stable_eq(j).x;
            y_eq = normU(x_eq);
            [~, idx] = min(abs(stable_eq(j).domain - x_eq));
            closest_point = stable_eq(j).domain(idx);
            width = abs(x_eq - closest_point);
            y_closest = normU(closest_point);
            depth = y_closest - y_eq;
            slope = depth / width;
            widths(j) = width;
            depths(j) = depth;
            slopes(j) = slope;
            eq{j} = [x_eq, y_eq];
        end
        all_slopes{k, i} = slopes;
        all_widths{k, i} = widths;
        all_depths{k, i} = depths;
        all_eq{k, i} = eq;
    end
end

%%
all_best_bins = zeros(1,n_subs);
for K = 1:n_subs
    for j = 1:40
        rmse(j) = all_mods{j, K}.xtra.reconstr.gofD2.rmse;
    end
    [~, idx_rmse] = sort(rmse);

    for j = 1:40
        ueff = all_mods{idx_rmse(j), K}.potential_eff;
        U = ueff.ueff(ueff.dom);
        if isreal(U)
            all_best_bins(K) = bins(idx_rmse(j));
            break
        end
    end
end
%% Analyses of some features

% 1. Gather derivatives from the potential plots to assess overall directionality
% 2. Depth of the eff_potential stable point/s to estimate the resilience
% 3. Width of the eff_potential stable point/s to estimate resilience
% 4. Slope of the eff_potential stable point/s
% 5. Number of stable points

median_crt = cellfun(@(v) median(v(:, 2)), sorted_tests_with_time);
median_srt = cellfun(@(v) median(v(:, 1)), sorted_tests_with_time);

%% no cheat
valid_age = [];
valid_pred_age = [];
valid_gender = [];
valid_depth = [];
valid_width = [];
valid_slope = [];
valid_dist0 = [];
valid_crt = [];
valid_srt = [];
valid_ranges = [];
valid_derivative = [];
valid_depth_sum = [];

for i = 1:n_subs
    if sorted_genders(i) == 999
        continue
    end
    if median_crt(i) == 0
        continue
    end
    if ~all_has_unstable(i)
        continue
    end
    % removing cause depth is ~50 massive outliers
    %     if i == 22 || i == 139
    %         continue
    %     end
    if ~isempty(all_depths{i})

        valid_age = [valid_age, sorted_ages(i)];
        valid_pred_age = [valid_pred_age, mean(sorted_preds{i})];
        valid_gender = [valid_gender, sorted_genders(i)];
        valid_crt = [valid_crt, median_crt(i)];
        valid_srt = [valid_srt, median_srt(i)];
        valid_derivative = [valid_derivative, derivative_patterns(i)];
        valid_ranges = [valid_ranges, all_dom_ranges(i) * all_SIGMA(i)];
        valid_depth_sum = [valid_depth_sum, sum(all_depths{i})];
        if length(all_depths{i}) > 1
            [~, idx] = max(all_depths{i});
            valid_depth = [valid_depth, all_depths{i}(idx)];
            valid_slope = [valid_slope, all_slopes{i}(idx) / all_SIGMA(i)];
            valid_width = [valid_width, all_widths{i}(idx) * all_SIGMA(i)];
            xy = all_eq{i}{idx};
            valid_dist0 = [valid_dist0, xy(1)];
        else
            valid_depth = [valid_depth, all_depths{i}];
            valid_slope = [valid_slope, all_slopes{i} / all_SIGMA(i)];
            valid_width = [valid_width, all_widths{i} * all_SIGMA(i)];
            xy = all_eq{i}{1};
            valid_dist0 = [valid_dist0, xy(1)];
        end
    end
end

tbl_no_cheat = array2table([valid_age; ...
    valid_pred_age; ...
    valid_age - valid_pred_age; ...
    valid_gender; ...
    valid_depth; ...
    valid_width; ...
    valid_slope; ...
    valid_crt; ...
    valid_srt; ...
    valid_dist0; ...
    valid_ranges; ...
    valid_depth_sum; ...
    valid_derivative]', 'VariableNames',{'age', 'pred_age','age_gap', 'gender', 'depth', 'width', 'slope', 'crt', 'srt', 'dist0','dom_range', 'depth_sum', 'derivative'});

%% cheat
valid_age = [];
valid_pred_age = [];
valid_gender = [];
valid_depth = [];
valid_width = [];
valid_slope = [];
valid_dist0 = [];
valid_crt = [];
valid_srt = [];
valid_ranges = [];
valid_derivative = [];
valid_depth_sum  = [];

for i = 1:n_subs
    if sorted_genders(i) == 999
        continue
    end
    if median_crt(i) == 0
        continue
    end
    %     if ~all_has_unstable(i)
    %         continue
    %     end
    % removing cause depth is ~50 massive outliers
    %     if i == 22 || i == 139
    %         continue
    %     end
    if ~isempty(all_depths{i})

        valid_age = [valid_age, sorted_ages(i)];
        valid_pred_age = [valid_pred_age, mean(sorted_preds{i})];
        valid_gender = [valid_gender, sorted_genders(i)];
        valid_crt = [valid_crt, median_crt(i)];
        valid_srt = [valid_srt, median_srt(i)];
        valid_derivative = [valid_derivative, derivative_patterns(i)];
        valid_ranges = [valid_ranges, all_dom_ranges(i) * all_SIGMA(i)];
        valid_depth_sum = [valid_depth_sum, sum(all_depths{i})];
        if length(all_depths{i}) > 1
            [~, idx] = max(all_depths{i});
            valid_depth = [valid_depth, all_depths{i}(idx)];
            valid_slope = [valid_slope, all_slopes{i}(idx) / all_SIGMA(i)];
            valid_width = [valid_width, all_widths{i}(idx) * all_SIGMA(i)];
            xy = all_eq{i}{idx};
            valid_dist0 = [valid_dist0, xy(1)];
        else
            valid_depth = [valid_depth, all_depths{i}];
            valid_slope = [valid_slope, all_slopes{i} / all_SIGMA(i)];
            valid_width = [valid_width, all_widths{i} * all_SIGMA(i)];
            xy = all_eq{i}{1};
            valid_dist0 = [valid_dist0, xy(1)];
        end
    end
end

tbl_cheat = array2table([valid_age; ...
    valid_pred_age; ...
    valid_age - valid_pred_age; ...
    valid_gender; ...
    valid_depth; ...
    valid_width; ...
    valid_slope; ...
    valid_crt; ...
    valid_srt; ...
    valid_dist0; ...
    valid_ranges; ...
    valid_depth_sum; ...
    valid_derivative]', 'VariableNames',{'age', 'pred_age','age_gap', 'gender', 'depth', 'width', 'slope', 'crt', 'srt', 'dist0', 'dom_range', 'depth_sum', 'derivative'});

%% Question 1 - Does the observation of resilient behavioural dynamics vary with age?
% no unstable points means there is only 1 "very" stable point, it seems
% this happens more often for older people

[H1,P1,CI1,STATS1] = ttest2(sorted_ages(all_has_unstable == 1), sorted_ages(all_has_unstable == 0));
mean_more_up = mean(sorted_ages(all_has_unstable == 1));
mean_no_up = mean(sorted_ages(all_has_unstable == 0));

%% Question 2 - Does the observation of two or more stable points vary with age?
% more stable points means more modalities, it seems this happens more for
% younger people
[H2,P2,CI2,STATS2] = ttest2(sorted_ages(n_stable_points == 1), sorted_ages(n_stable_points > 1));
mean_one_sp = mean(sorted_ages(n_stable_points == 1));
mean_more_sp = mean(sorted_ages(n_stable_points > 1));

%% Question 3 - Do the measures of the landscape and resiliance vary with age?

mdl_3_no_cheat = fitlm(tbl_no_cheat, 'age ~ gender + derivative + slope + width + depth', 'RobustOpts', 'on');
mdl_3_cheat = fitlm(tbl_cheat, 'age ~ gender + derivative + slope + width + depth', 'RobustOpts', 'on');

%% Question 4 - followup does the derivative vary with age

mdl_4_cheat = fitlm(tbl_cheat, 'age ~ gender + derivative', 'RobustOpts', 'on');
mdl_x_cheat = fitlm(tbl_cheat, 'age_gap ~ gender + derivative', 'RobustOpts', 'on');

%% Question 5 - Does the task performance (crt) depend on the behavioral dyanmics?

mdl_5_no_cheat = fitlm(tbl_no_cheat, 'crt ~ age + gender + derivative + slope + width + depth', 'RobustOpts', 'on');
mdl_5_cheat = fitlm(tbl_cheat, 'crt ~ age + gender + derivative', 'RobustOpts', 'on');

%% Question 6 - Does the task performance (srt) depend on the behavioral dyanmics?

mdl_6_no_cheat = fitlm(tbl_no_cheat, 'srt ~ age + gender + derivative + slope + width + depth', 'RobustOpts', 'on');
mdl_6_cheat = fitlm(tbl_cheat, 'srt ~ age + gender + derivative', 'RobustOpts', 'on');


%% UTILS
%% plot after
[~, idx_der] = sort(derivative_patterns);
K = idx_der(end-7);
x_ = double(all_preds{idxs(K)});
SIGMA = all_SIGMA(K);
MU = all_MU(K);
ueff = all_mods{1, K}.potential_eff;
U = ueff.ueff(ueff.dom);
U = U - min(U) + 0.1 * (max(U) - min(U));
new_dom = ueff.dom * SIGMA + MU;

subplot(3,1,2)
plot(new_dom, U)
hold on
yl = ylim;
xlim([min(new_dom), max(new_dom)])
% plot([sorted_ages(K), sorted_ages(K)], yl, 'r')
title(sprintf('%f   %f', sorted_ages(K), derivative_patterns(i)))
% plot([mean(x_), mean(x_)], yl, 'b')
eq_points = all_eq{K};
for i = 1:length(eq_points)
    x_eq = eq_points{i}(1) * SIGMA + MU;
    y_eq = eq_points{i}(2);
    new_mu = all_widths{K}(i) * SIGMA;
    scatter(x_eq, y_eq, 'c')
    plot(sort([x_eq - new_mu,  x_eq]), [y_eq, y_eq], 'g')
    plot([x_eq - new_mu,x_eq - new_mu], sort([y_eq, y_eq + all_depths{K}(i)]), 'b')
end
xl = xlim;
subplot(3,1,1)
plot(all_preds{idxs(K)}, 1:length(all_preds{idxs(K)}))
xlim(xl)
set(gca, 'YDir','reverse')
subplot(3,1,3)
u = all_mods{1, K}.potential;
U = u.u(u.dom);
U = U - min(U) + 0.1 * (max(U) - min(U));
new_dom = u.dom * SIGMA + MU;
plot(new_dom, U)
xlim([min(new_dom), max(new_dom)])
