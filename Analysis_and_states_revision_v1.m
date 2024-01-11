load('all_mean_exit_times_28_12_2023.mat')
load('analysis_output_landscapes_15_12_2023.mat')


%% Eliminations 

% Elimantate individuals who have no stable points 
idx_nostable = (tbl.n_stable_points == 0);
tbl(idx_nostable,:) = [];
all_exit_times(idx_nostable) = []; 
% Eliminate individuals with too many wells (as in https://doi.org/10.1371/journal.pone.0041010 Such high number of detected states meant that, in principle, the data were on the edge of having no clear potential) 
idx_nustable = (tbl.n_stable_points > 5);
tbl(idx_nustable,:) = [];
all_exit_times(idx_nustable) = []; 
% Mark people who did not report gender with NaN instead of 999
tbl.gender(tbl.gender>2) = NaN;%

%% 

n_unstable_points = tbl.n_unstable_points;

%% Description: 
figure; title ('Number of unstable points'); hold on; 
ecdf(tbl.n_unstable_points) % [***2C***]

%% Linking overall propensity of the system 
tmp_derivative =  tbl.derivative; 
tmp_derivative(tmp_derivative>0.2) = NaN;   % remove outliers 
%Mdl.derivative = fitlm([tbl.age], tmp_derivative, 'VarNames', {'Age' ,'Derivative'}, 'RobustOpts', 'on')

tbl.age_gap(tbl.age_gap <-40) = NaN; 
%[Figure with age]
Mdl.derivative = fitlm([-tbl.age_gap tbl.age],tmp_derivative,'VarNames' , {'Age_gap' ,'Age', 'Derivative'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***2D*** focused on Age on X axis]

handles2d = plotAdjustedResponse(Mdl.derivative, 'Age'); 
adjusted2d = prepare_handle_for_python(handles2d);
adjusted2d.pval = Mdl.derivative.Coefficients{'Age', 'pValue'};

save('for_figure_2', 'n_unstable_points', 'adjusted2d')
%% Case 1: One unstable point/ Case 2: More than one unstable point 
% go through all subjects and find the young_to_old
% and old_to_young [y2o & o2y] exit times. In case of one unstable point
% 'AR' > old to young : 'RA' > young to old 

for i = 1:size(tbl,1)

    tmp_e = all_exit_times{i};
    tmp_exit_y2o = NaN;
    tmp_exit_o2y = NaN;
    tmp_l = 0; tmp_m = 0;
    for t = 1:length(tmp_e)
        if strcmp(tmp_e{t}.BC, 'AR')
            tmp_l = tmp_l+1;
            tmp_exit_o2y(tmp_l) = tmp_e{t}.WT;
        end
        if strcmp(tmp_e{t}.BC, 'RA')
            tmp_m = tmp_m+1;
            tmp_exit_y2o(tmp_m) = tmp_e{t}.WT;
        end

    end
    exit_y2o(i) = mean(tmp_exit_y2o); % there should be only one in any case
    exit_o2y(i) = mean(tmp_exit_o2y);  % there should be only one in any case
    %end
    clear tmp*
end
exit_y2o(exit_y2o==0) = deal(NaN);
exit_o2y(exit_o2y==0) = deal(NaN); 

% young to old [Figure with Age, insert ONE]
Mdl.y2o = fitlm([-tbl.age_gap tbl.age],log10(exit_y2o),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***3B***]
Mdl.y2o_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points==1) tbl.age(tbl.n_unstable_points==1)],log10(exit_y2o(tbl.n_unstable_points==1)),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***Suppl Fig 7 ***]

handles3b = plotAdjustedResponse(Mdl.y2o, 'Age'); 
adjusted3b = prepare_handle_for_python(handles3b);
adjusted3b.pval = Mdl.y2o.Coefficients{'Age', 'pValue'};

handles7b = plotAdjustedResponse(Mdl.y2o_ONE, 'Age'); 
adjusted7b = prepare_handle_for_python(handles7b);
adjusted7b.pval = Mdl.y2o_ONE.Coefficients{'Age', 'pValue'};


% old to young [Figure with Age, insert ONE]
Mdl.o2y = fitlm([-tbl.age_gap tbl.age],log10(exit_o2y),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***3C***]
Mdl.o2y_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points==1) tbl.age(tbl.n_unstable_points==1)],log10(exit_o2y(tbl.n_unstable_points==1)),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***Suppl Fig 7 ***]

handles3c = plotAdjustedResponse(Mdl.o2y, 'Age'); 
adjusted3c = prepare_handle_for_python(handles3c);
adjusted3c.pval = Mdl.o2y.Coefficients{'Age', 'pValue'};

handles7c = plotAdjustedResponse(Mdl.o2y_ONE, 'Age'); 
adjusted7c = prepare_handle_for_python(handles7c);
adjusted7c.pval = Mdl.o2y_ONE.Coefficients{'Age', 'pValue'};

save('for_figure_3_bc.mat', 'adjusted7c', 'adjusted7b', 'adjusted3c', 'adjusted3b')

% gap between young to old and old to young [set up as o2y - y2o, so higher
% the value the worse] [Text results only]
exit_gap = exit_o2y-exit_y2o; 
Mdl.exitgap = fitlm([-tbl.age_gap tbl.age],exit_gap,'VarNames' , {'Age_gap' ,'Age', 'Exit gap'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse)
Mdl.exitgap_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points==1) tbl.age(tbl.n_unstable_points==1)],exit_gap(tbl.n_unstable_points==1),'VarNames' , {'Age_gap' ,'Age', 'Exit gap'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse)

%% Address where are the alternative stable states relative to the most stable point? 
idx_left_tmp = tbl.depth_left>0;
idx_right_tmp = tbl.depth_right>0;
idx_both = [idx_right_tmp&idx_left_tmp];

idx_left = [idx_left_tmp&~idx_both]; 
idx_right = [idx_right_tmp&~idx_both];

figure; title('Is the younger alternative available irrespective of age?'); %[****Suppl Fig. 6***** ]
ecdf(tbl.age(idx_left)); hold on
ecdf(tbl.age(idx_right))
legend({'Younger alt', 'Older alt'})
% people who have an younger alternative are typically older than the people who have
% an older alternative

age = tbl.age;

save('for_figure_2', 'n_unstable_points', 'adjusted2d', 'age', 'idx_right', 'idx_left')

% clearvars -except Mdl tbl


%% introduce mean exit time
% 127 looks good
% 74, 76
% 65, 68, 88 also good

idx = 127;
x = sorted_preds{idx};
t = sorted_pred_times{idx};
MU = all_MU(idx);
SIGMA = all_SIGMA(idx);
subplot(2,1,1)
mod = all_mods{idx};
mod.equilibria = mod.find_equilibria('effective');
mod.plot('potential_eff')

dom = mod.potential_eff.dom;
ueff = mod.potential_eff.ueff;
U = ueff(mod.potential_eff.dom);
normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));
y = normU(dom);

subplot(2,1,2)
mean_exit_all = mod.mean_exit('all', mod.pdf);
mod.plot('mean_exit')

met_T_1 = mean_exit_all{1}.T;
met_T_2 = mean_exit_all{2}.T;
met_WT_1 = mean_exit_all{1}.WT;
met_WT_2 = mean_exit_all{2}.WT;
met_dom_2 = mean_exit_all{2}.domain;
met_dom_1 = mean_exit_all{1}.domain;

save('for_figures_3_4', 'x', 't', ...
    'y', 'dom', 'MU', 'SIGMA', ...
    'met_T_2', 'met_T_1', 'met_WT_2', 'met_WT_1', ...
    'met_dom_1', 'met_dom_2')

%% prepare handle for python

function out_handle = prepare_handle_for_python(in_handle)

    adjusted_data = zeros(2, size(in_handle(1,1).XData, 2));
    adjusted_fit = zeros(2, size(in_handle(2,1).XData, 2));
    adjusted_data(1, :) = in_handle(1,1).XData;
    adjusted_data(2, :) = in_handle(1,1).YData;
    adjusted_fit(1, :) = in_handle(2,1).XData;
    adjusted_fit(2, :) = in_handle(2,1).YData;
    
    out_handle.adjusted_fit = adjusted_fit;
    out_handle.adjusted_data = adjusted_data;
end