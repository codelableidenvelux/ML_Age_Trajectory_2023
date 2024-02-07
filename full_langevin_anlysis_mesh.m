%% load healthy age data 20d
load('all_age_pred_21d_30_11_2023.mat')

[sorted_ages, idxs] = sort(cellfun(@(x) x(end), all_ages));
sorted_genders = all_genders(idxs);
sorted_preds = all_preds(idxs);
sorted_pred_times = all_preds_times(idxs);

sorted_pids = all_pids(idxs);

%% run analysis

% prepare variables
n_subs = length(all_ages);
derivative_patterns = zeros(1, n_subs);
median_predicted_ages = zeros(1, n_subs);
gap_real_age_median_predicted_age = zeros(1, n_subs);
n_equilibria = zeros(1, n_subs);
n_stable_points = zeros(1, n_subs);
n_unstable_points = zeros(1, n_subs);
all_MU = zeros(1, n_subs);
all_SIGMA = zeros(1, n_subs);
all_slopes = cell(1, n_subs);
all_depths = cell(1, n_subs);
all_widths = cell(1, n_subs);
all_dists = zeros(1, n_subs);
all_dom_ranges = zeros(1, n_subs);
all_eq = cell(1, n_subs);
all_mods = cell(1, n_subs);
all_has_unstable = zeros(1, n_subs);
all_left_right_combo = zeros(4, n_subs);
all_left_right_depth = zeros(2, n_subs);
all_left_right_width = zeros(2, n_subs);
all_exit_times = cell(1, n_subs);
exit_times = zeros(1, n_subs);
all_stable_points = cell(1, n_subs);


% We have 4 cases
% 1 - 0 unstable points mean we have left edge slope and right edge slope
% 2 - at least 2 unstable points one right and one left of the deepest stable point
% 3 - exactly 1 unstable point on the left -> left unstable slope and right edge slope
% 4 - extacly 1 unstable point on the right -> left edge slope, right unstable slope

for i = 52:n_subs
    fprintf("%d\n", i);
    % get trajectory for subject
    x_ = double(sorted_preds{i});

    % calculate median predicted age over the trajectory
    median_predicted_ages(i) = median(x_, 'omitnan');

    % calculate the gap between the median prediction and the real age
    gap_real_age_median_predicted_age(i) = median(x_, 'omitnan') - sorted_ages(i);

    % z-scoring the trajectories - allows us to come back to the age scale
    % later
    MU = mean(x_, 'omitnan');
    SIGMA = std(x_, 'omitnan');
    x_ = (x_ - MU) / SIGMA;
    all_MU(i) = MU;
    all_SIGMA(i) = SIGMA;

    avec = linspace(min(x_), max(x_), length(x_) / 5);
    bw = 0.3 * std(x_, 'omitnan');
    DT = 1;
    nx = length(diff(x_));
    % Calculate the Langevin reconstruction by means of the mesh method 
    % as opposed to the bin method in Arani et al.
    results_M = LangevinReconst_MESH(x_,diff(x_),nx,DT,bw,length(avec),avec);

    % get the model from the results
    mod = langevin_eq(results_M);
    all_mods{i} = mod;

    % use functionality to calculate potential eff
    ueff = mod.potential_eff.ueff;
    U = ueff(mod.potential_eff.dom);
    normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));

    % sometimes we get complex numbers, model has failed, we skip sample
    if ~isreal(U)
        continue
    end

    % first order approximation of the average derivative of the DT landscape
    derivative_x =  diff(U);
    derivative_patterns(i) = mean(derivative_x);

    all_dom_ranges(i) = mod.potential_eff.dom(end) - mod.potential_eff.dom(1);

    % extract depth - width - slope
    % calculate equilibria of the effective potential landascape
    mod.equilibria = mod.find_equilibria('effective');

    n_equilibria(i) = length(mod.equilibria);

    % extract stable and unstable points
    unstable_eq = mod.equilibria([mod.equilibria.stable] == 0);
    stable_eq = mod.equilibria([mod.equilibria.stable] == 1);

    all_stable_points{i} = stable_eq;

    n_stable_points(i) = length(stable_eq);
    n_unstable_points(i) = length(unstable_eq);

    % wheter or not an unstable point is present
    all_has_unstable(i) = ~isempty(unstable_eq);

    % pre-allocation
    slopes = zeros(1, length(stable_eq));
    eq = cell(1, length(stable_eq));

    % we do not consider cases (1) when no stable points are present
    if isempty(stable_eq)
        continue
    end

    % no unstable points: case 1 with only 1 stable point. We clculate
    % slope of left and right edges.
    if (~isempty(stable_eq)) && (isempty(unstable_eq))

%         mean_exit_all = mod.mean_exit('all', mod.pdf);
%         all_exit_times{i} = mean_exit_all;

        x_eq = stable_eq(1).x;
        y_eq = normU(x_eq);

        x_left = mod.potential_eff.dom(2);
        x_right = mod.potential_eff.dom(end-1);

        y_left = normU(x_left);
        y_right = normU(x_right);

        slope_left = calc_slope(x_eq, y_eq, x_left, y_left);
        slope_right = calc_slope(x_eq, y_eq, x_right, y_right);

        all_left_right_combo([1, 4], i) = [slope_left, slope_right];
        all_eq{i} = {[x_eq, y_eq]};

    elseif length(unstable_eq) == 1 % here is case 3 or 4

        % calculate mean exit time for all stable points
        mean_exit_all = mod.mean_exit('all', mod.pdf);
        all_exit_times{i} = mean_exit_all;
        % after this step one can plot the mean exit time
        % mod.plot('mean_exit')
        % mean exit time is calculated for each of the stable points
        % the calculation of each mean exit time depend on what is on the 
        % left and what is on the right of the stable point, i.e. an
        % unstable point or an edge. This can be seen from 
        % mean_exit_all{i}.BC where 'R' means repulsive generally an edge
        % and 'A' means attractive generally an unstable point.
        % When you have 'RA' or 'AR' generally the mean exit time is
        % simpler to calculate. When you have 'AA' you are in between 2
        % unstable points, this means that the mean exit time is a little
        % more complicated, thus the calculation includes attraction from
        % both left and right side which you can still see in the plot in
        % blue and magenta. At the end the final value of mean exit time is
        % always contained in mean_exit_all{i}.WT and is expressed in days.
        % to calculate



        % find most stable point: highest mean exit time. The most stable
        % point can still be in case 3 or 4, meaning having the unstable 
        % point on the right and on the left.
        [~, long_exit] = max(cellfun(@(x) x.WT, mean_exit_all));
        exit_times(i) = mean_exit_all{long_exit}.WT;

        x_eq = stable_eq(long_exit).x;
        y_eq = normU(x_eq);

        x_un = unstable_eq(1).x;
        y_un = normU(x_un);

        x_left = mod.potential_eff.dom(2);
        x_right = mod.potential_eff.dom(end-1);

        y_left = normU(x_left);
        y_right = normU(x_right);

        [slope_left, width_left, depth_left] = calc_slope(x_eq, y_eq, x_left, y_left);
        [slope_right, width_right, depth_right] = calc_slope(x_eq, y_eq, x_right, y_right);
        [slope, width, depth] = calc_slope(x_eq, y_eq, x_un, y_un);

        if slope > 0 % unstable point on the right - 3
            all_left_right_combo([1, 3], i) = [slope_left, slope];
            all_left_right_width(2, i) = width;
            all_left_right_depth(2, i) = depth;
        else % unstable point on the right - 4
            all_left_right_combo([2, 4], i) = [slope, slope_right];
            all_left_right_width(1, i) = width;
            all_left_right_depth(1, i) = depth;
        end

        all_eq{i} = {[x_eq, y_eq]};

    else  % multiple stable and unstable points

        % calculate mean exit time for all stable points
        mean_exit_all = mod.mean_exit('all', mod.pdf);
        all_exit_times{i} = mean_exit_all;

        % find most stable point (max mean exit time)
        [~, long_exit] = max(cellfun(@(x) x.WT, mean_exit_all));
        exit_times(i) = mean_exit_all{long_exit}.WT;

        x_eq = stable_eq(long_exit).x;
        y_eq = normU(x_eq);

        x_un = unstable_eq(1).x;
        y_un = normU(x_un);

        % here I need to check if the domain on the right or on the left of 
        % the unstable point is Inf, if it is it means in that direction I 
        % have an edge, otherwise I have a stable point  

        % left domain
        if abs(stable_eq(long_exit).domain(1)) == Inf
            x_left = mod.potential_eff.dom(2);
            y_left = normU(x_left);
            [slope_left, width_left, depth_left] = calc_slope(x_eq, y_eq, x_left, y_left);
        else
            x_un = stable_eq(long_exit).domain(1);
            y_un = normU(x_un);
            [slope_left, width_left, depth_left] = calc_slope(x_eq, y_eq, x_un, y_un);
        end

        % right domain
        if abs(stable_eq(long_exit).domain(2)) == Inf
            x_right = mod.potential_eff.dom(end-1);
            y_right = normU(x_right);
            [slope_right, width_right, depth_right] = calc_slope(x_eq, y_eq, x_right, y_right);
        else
            x_un = stable_eq(long_exit).domain(2);
            y_un = normU(x_un);
            [slope_right, width_right, depth_right] = calc_slope(x_eq, y_eq, x_un, y_un);
        end

        if (abs(stable_eq(long_exit).domain(1)) ~= Inf) && (abs(stable_eq(long_exit).domain(2)) ~= Inf)
            all_left_right_combo([2, 3], i) = [slope_left, slope_right];
            all_left_right_depth([1, 2], i) = [depth_left, depth_right];
            all_left_right_width([1, 2], i) = [width_left, width_right];
        elseif (abs(stable_eq(long_exit).domain(1)) ~= Inf)
            all_left_right_combo([2, 4], i) = [slope_left, slope_right];
            all_left_right_depth(1, i) = depth_left;
            all_left_right_width(1, i) = width_left;
        else
            all_left_right_combo([1, 3], i) = [slope_left, slope_right];
            all_left_right_depth(2, i) = depth_right;
            all_left_right_width(2, i) = width_right;
        end
    end
    if ~isreal(all_left_right_combo(:, i))
        fprintf('mmm')
    end
end

% given the assumptions on what comes first some value might be negative 
% but it should always treated as positive
all_left_right_depth = abs(all_left_right_depth);
all_left_right_width = abs(all_left_right_width);

%% assemble data

% mean_predicted_age = cellfun(@(x) mean(x), sorted_preds);
% median_crt = cellfun(@(v) median(v(:, 2)), sorted_tests_with_time);
% median_srt = cellfun(@(v) median(v(:, 1)), sorted_tests_with_time);


tbl = array2table([sorted_ages; ...
    median_predicted_ages; ...
    sorted_ages - median_predicted_ages; ...
    cell2mat(sorted_genders); ...
    exit_times; ...
    all_left_right_depth; ...
    all_left_right_width; ...
    all_left_right_combo; ...
%     median_crt; ...
%     median_srt; ...
    all_dom_ranges; ...
    all_MU ; ...
    all_SIGMA; ...
    n_unstable_points; ...
    n_stable_points; ...
    derivative_patterns]', ...
    'VariableNames',{ ...
    'age', ...
    'pred_age', ...
    'age_gap', ...
    'gender', ...
    'mean_exit_times', ...
    'depth_left', 'depth_right', ...
    'width_left', 'width_right', ...
    'slope_edge_left', 'slope_unstable_left', 'slope_unstable_right', 'slope_edge_right', ...
    'dom_range', 'mu', 'sigma', ...
    'n_unstable_points', 'n_stable_points', ...
    'derivative'});


%% double check -> just visual checking that the above is correct 

for i = 1:n_subs
    close all
    all_mods{i}.plot('potential_eff')
    ueff = all_mods{i}.potential_eff.ueff;
    U = ueff(all_mods{i}.potential_eff.dom);
    normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));

    stable_eq = all_mods{i}.equilibria([all_mods{i}.equilibria.stable] == 1);
    if isempty(stable_eq)
        continue
    end

    for j = 1:4
        mean_exit_all = mod.mean_exit('all', mod.pdf);
   
        % find most stable point (max mean exit time)
        [~, long_exit] = max(cellfun(@(x) x.WT, mean_exit_all));
        s = all_left_right_combo(j, i);
        x1 = stable_eq(long_exit).x;
        y1 = normU(x1);
        if s > 0
            x2 = all_mods{i}.potential_eff.dom(end-1);
            y2 = y1 + s * (x2 - x1);
            plot([x1, x2], [y1, y2], '--r', 'LineWidth', 3)
        elseif s < 0
            x2 = all_mods{i}.potential_eff.dom(2);
            y2 = y1 + s * (x2 - x1);
            plot([x2, x1], [y2, y1], '--r', 'LineWidth', 3)
        end

    end

    % left

    [~, deep_idx] = min(cellfun(normU, {stable_eq.x}));
    s = all_left_right_combo(j, i);
    x1 = stable_eq(deep_idx).x;
    y1 = normU(x1);

    w = all_left_right_width(1, i);
    d = all_left_right_depth(1, i);

    if w ~= 0
        plot([x1 - w, x1], [y1, y1], '--g', 'LineWidth', 3)
        plot([x1 - w, x1 - w], [y1, y1 + d], '--y', 'LineWidth', 3)
    end

    % right
    w = all_left_right_width(2, i);
    d = all_left_right_depth(2, i);

    if w ~= 0
        plot([x1 + w, x1], [y1, y1], '--g', 'LineWidth', 3)
        plot([x1 + w, x1 + w], [y1, y1 + d], '--y', 'LineWidth', 3)
    end


    title(i)
    pause;

end


%% plot after -> this is for finding examples for figure 2
[~, idx_der] = sort(derivative_patterns);

% this_pid = '138ebb8bf527057147e2b63aa6b026535b2028eb';
% this_pid = '138e9cefa56279c54afda02470df4da23dcd28eb';
% this_pid = '138eb833b8fd25944a228eb31a67909eace928eb';
% this_pid = '138ed7413bc320804ce9923caf034bc24cb428eb'; 
% this_pid = '138eb8d6ccaacc9742dc836c430cbbc0dc7828eb';
% this_pid = '138ea57b58d081e8485f873d1e670454a88928eb';
this_pid = '138ebfe89d81730c41db82f2accec7286e8928eb';


for i = 1:265
    if (sorted_pids(i,:) == this_pid)
    K = i;
    end
end

x_ = double(all_preds{idxs(K)});
SIGMA = all_SIGMA(K);
MU = all_MU(K);
real_age = sorted_ages(K);
ueff = all_mods{1, K}.potential_eff;
U = ueff.ueff(ueff.dom);
U = U - min(U) + 0.1 * (max(U) - min(U));
new_dom = ueff.dom * SIGMA + MU;

subplot(2,1,2)
plot(new_dom, U)
hold on
yl = ylim;
xlim([min(new_dom), max(new_dom)])
% plot([sorted_ages(K), sorted_ages(K)], yl, 'r')
title(sprintf('%f   %f', sorted_ages(K), derivative_patterns(i)))
plot([mean(x_, 'omitnan'), mean(x_,  'omitnan')], yl, 'b--')
plot([real_age, real_age], yl, 'r')
% eq_points = all_eq{K};
% depths = all_depths{K};
% widths = all_widths{K};
% slopes = all_slopes{K};
% for i = 1:length(eq_points)
%     x_eq = eq_points{i}(1) * SIGMA + MU;
%     y_eq = eq_points{i}(2);
%     new_mu = all_widths{K}(i) * SIGMA;
%     scatter(x_eq, y_eq, 'c')
%     xxlim = xlim;
%     if (x_eq - new_mu > xxlim(1))
%         plot(sort([x_eq - new_mu,  x_eq]), [y_eq, y_eq], 'g')
%         plot([x_eq - new_mu,x_eq - new_mu], sort([y_eq, y_eq + all_depths{K}(i)]), 'b')
%     else
%         plot(sort([x_eq + new_mu,  x_eq]), [y_eq, y_eq], 'g')
%         plot([x_eq + new_mu,x_eq + new_mu], sort([y_eq, y_eq + all_depths{K}(i)]), 'b')
%     end
% end
xl = xlim;
subplot(2,1,1)
plot(all_preds{idxs(K)}, 1:length(all_preds{idxs(K)}))
xlim(xl)
set(gca, 'YDir','reverse')


i = K;
figure(2)
all_mods{i}.plot('potential_eff')
ueff = all_mods{i}.potential_eff.ueff;
U = ueff(all_mods{i}.potential_eff.dom);
normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));

stable_eq = all_mods{i}.equilibria([all_mods{i}.equilibria.stable] == 1);


for j = 1:4
    [~, deep_idx] = min(cellfun(normU, {stable_eq.x}));
    s = all_left_right_combo(j, i);
    x1 = stable_eq(deep_idx).x;
    y1 = normU(x1);
    if s > 0
        x2 = all_mods{i}.potential_eff.dom(end-1);
        y2 = y1 + s * (x2 - x1);
        plot([x1, x2], [y1, y2], '--r', 'LineWidth', 3)
    elseif s < 0
        x2 = all_mods{i}.potential_eff.dom(2);
        y2 = y1 + s * (x2 - x1);
        plot([x2, x1], [y2, y1], '--r', 'LineWidth', 3)
    end

end

% left

[~, deep_idx] = min(cellfun(normU, {stable_eq.x}));
s = all_left_right_combo(j, i);
x1 = stable_eq(deep_idx).x;
y1 = normU(x1);

eq_points = [x1, y1];

w = all_left_right_width(1, i);
d = all_left_right_depth(1, i);

if w ~= 0
    plot([x1 - w, x1], [y1, y1], '--g', 'LineWidth', 3)
    plot([x1 - w, x1 - w], [y1, y1 + d], '--y', 'LineWidth', 3)
end

% right
w = all_left_right_width(2, i);
d = all_left_right_depth(2, i);

if w ~= 0
    plot([x1 + w, x1], [y1, y1], '--g', 'LineWidth', 3)
    plot([x1 + w, x1 + w], [y1, y1 + d], '--y', 'LineWidth', 3)
end


title(i)

U = normU(all_mods{i}.potential_eff.dom);
depths = all_left_right_depth(:, i);
widths = all_left_right_width(:, i);

save(sprintf('example_idx_%s.mat', this_pid), 'real_age', 'MU', 'widths', 'depths', 'SIGMA', 'U', 'eq_points', 'new_dom', 'x_')


%% utility function to calculate slope
function [slope, width, depth] = calc_slope(x1, y1, x2, y2)
    width = x1 - x2;
    depth = y1 - y2;
    slope = depth / width;
end


