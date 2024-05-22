load('all_mean_exit_times_28_12_2023.mat')
load('analysis_output_landscapes_15_12_2023.mat')
x = readtable('median_taps_vs_chrono_age_plus_tregectories.csv'); 

%% Add duration to table 
for i = 1:size(all_preds,2)
tbl.dur(i) = size(all_preds{1,i},2);
end

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

%% Relationship between usage and age 

usage = fitlm(x.age(ismember(x.trajectory, 'yes')), sqrt(x.median_taps_per_day(ismember(x.trajectory, 'yes'))), 'RobustOpts', 'on'); 

%% Description: 
figure; title ('Number of unstable points'); hold on; 
ecdf(tbl.n_unstable_points) % [***2C***]

% number of people with at least one tipping point 



%% Linking overall propensity of the system 
tmp_derivative =  tbl.derivative; 
tmp_derivative(tmp_derivative>0.2) = NaN;   % remove outliers 
%Mdl.derivative = fitlm([tbl.age], tmp_derivative, 'VarNames', {'Age' ,'Derivative'}, 'RobustOpts', 'on')

tbl.age_gap(tbl.age_gap <-40) = NaN; 
%[Figure with age]
Mdl.derivative = fitlm([-tbl.age_gap tbl.age],tmp_derivative,'VarNames' , {'Age_gap' ,'Age', 'Derivative'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***2D*** focused on Age on X axis]


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
Mdl.y2o_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2) tbl.age(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2)],log10(exit_y2o(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2)),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on');

% old to young [Figure with Age, insert ONE]
Mdl.o2y = fitlm([-tbl.age_gap tbl.age],log10(exit_o2y),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***3C***]
Mdl.o2y_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2) tbl.age(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2)],log10(exit_o2y(tbl.n_unstable_points == 1 & tbl.n_stable_points == 2)),'VarNames' , {'Age_gap' ,'Age', 'Exit time'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***Suppl Fig 7 ***]

% gap between young to old and old to young [set up as o2y - y2o, so higher
% the value the worse] [Text results only]
exit_gap = exit_o2y-exit_y2o; 
Mdl.exitgap = fitlm([-tbl.age_gap tbl.age],exit_gap,'VarNames' , {'Age_gap' ,'Age', 'Exit gap'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse)
Mdl.exitgap_ONE = fitlm([-tbl.age_gap(tbl.n_unstable_points==1) tbl.age(tbl.n_unstable_points==1)],exit_gap(tbl.n_unstable_points==1),'VarNames' , {'Age_gap' ,'Age', 'Exit gap'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse)

%% Descriptive stats on exit time: 
allexit = nanmax(exit_y2o, exit_o2y);

Mdl.allexitmax = fitlm([-tbl.age_gap tbl.age],log10(allexit),'VarNames' , {'Age_gap' ,'Age', 'Derivative'},'RobustOpts', 'on'); % add sign to make pred - real (i.e., larger the gap the worse) [***2D*** focused on Age on X axis]

medianexit = nanmedian(allexit); 
iqrexit = iqr(allexit);

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

clearvars -except Mdl tbl