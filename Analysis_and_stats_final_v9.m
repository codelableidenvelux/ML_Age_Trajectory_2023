
%% load the data on aging 
load('organized_trajectory_analysys_output_20d_v9.mat');
load('all_pred_seen_and_unseen_with_crt_20d_25_04_2023.mat');

%% elimination and correction 
for i = 1:length(all_preds_times)
tbl.dur(i) = length(all_preds_times{1,i}); 
end
tbl(tbl.n_stable_points == 0,:) =[]; 
tbl.slope_edge_left(tbl.slope_edge_left == 0) = NaN; 
tbl.slope_unstable_left(tbl.slope_unstable_left == 0) = NaN; 
tbl.gender(tbl.gender>2) = NaN;

n_unstable_points = tbl.n_unstable_points;
age = tbl.age;

%% ECDF: Distribution of number of unstable points [Fig. 3B]
figure; title ('Number of unstable points'); hold on; 
ecdf(tbl.n_unstable_points)

%% VIOLIN PLOT: Q1: Is the median age of people with no unstable point vs. at least one unstable point any different? [Fig. 3C]
[H2,P2,~,STATS2] = ttest2(tbl.age(tbl.n_unstable_points == 0), tbl.age(tbl.n_unstable_points >= 1));
figure; title('Age, no unstable point vs. >= 1 unstable point'); hold on; 
scatter(zeros(length(tbl.age(tbl.n_unstable_points == 0))), tbl.age(tbl.n_unstable_points == 0))
hold on
scatter(ones(length(tbl.age(tbl.n_unstable_points >= 1))), tbl.age(tbl.n_unstable_points >= 1))
Val2_0 = tbl.age(tbl.n_unstable_points == 0);
Val2_1 = tbl.age(tbl.n_unstable_points >= 1);



%% SCATTER PLOT 1: How does the general tendency (derivative) alter as a function of age, gender 
tmp_derivative =  tbl.derivative; 
tmp_derivative(tmp_derivative<-0.2) = NaN;   % remove outliers 
Mdl.derivative = fitlm([tbl.age], tmp_derivative, 'VarNames', {'Age' ,'Derivative'}, 'RobustOpts', 'on');
%figure;  hold on; % [Fig. 4A]
%plotAdjustedResponse(Mdl.derivative, 'Age'); title ('Derivatives vs. age');
handles4d = plotAdjustedResponse(Mdl.derivative, 'Age'); title ('Derivatives vs. age');
adjusted4d = prepare_handle_for_python(handles4d);
adjusted4d.pval = Mdl.derivative.Coefficients{'Age', 'pValue'};


%% How hard is it get younger? (unstable point on the left side of the deepest valley) 

% from perspective of the slope 
tmp_slopeL = [tbl.slope_unstable_left];
tmp_slopeL(tmp_slopeL == 0) = NaN; 
Mdl.slopeL = fitlm([tbl.age], sqrt(-tmp_slopeL), 'VarNames', {'Age' , 'SlopeL'}, 'RobustOpts', 'on');
%plotAdjustedResponse(Mdl.slopeL, 'Age'); title ('Slope to younger vs. age');
handles4f = plotAdjustedResponse(Mdl.slopeL, 'Age'); title ('Slope to younger vs. age');
adjusted4f = prepare_handle_for_python(handles4f);
adjusted4f.pval = Mdl.slopeL.Coefficients{'Age', 'pValue'};

% from the perspective of the depth 
tmp_depthL = [tbl.depth_left];
tmp_depthL(tmp_depthL == 0) = NaN; 
Mdl.depthL = fitlm([tbl.age], sqrt(tmp_depthL), 'VarNames', {'Age' , 'DepthL'}, 'RobustOpts', 'on');
handles4h = plotAdjustedResponse(Mdl.depthL, 'Age'); title ('Depth to younger vs. age');
adjusted4h = prepare_handle_for_python(handles4h);
adjusted4h.pval = Mdl.depthL.Coefficients{'Age', 'pValue'};

% from the perspective of the width 
tmp_widthL = [tbl.width_left];
tmp_widthL(tmp_widthL == 0) = NaN; 
Mdl.widthL = fitlm([tbl.age], tmp_widthL, 'VarNames', {'Age' , 'WidthL'}, 'RobustOpts', 'on');
handles4j = plotAdjustedResponse(Mdl.widthL, 'Age'); title ('width to younger vs. age');
adjusted4j = prepare_handle_for_python(handles4j);
adjusted4j.pval = Mdl.widthL.Coefficients{'Age', 'pValue'};

%% How hard is it to get older? (unstable point on the right) 
tmp_slopeR = [tbl.slope_unstable_right];
tmp_slopeR(tmp_slopeR == 0) = NaN; 
Mdl.slopeR = fitlm([tbl.age], sqrt(tmp_slopeR), 'VarNames', {'Age' , 'SlopeR'}, 'RobustOpts', 'on');
%plotAdjustedResponse(Mdl.slopeLU, 'Age'); title ('Slope to younger vs. age');
handles4e = plotAdjustedResponse(Mdl.slopeR, 'Age'); title ('Slope to younger vs. age');
adjusted4e = prepare_handle_for_python(handles4e);
adjusted4e.pval = Mdl.slopeR.Coefficients{'Age', 'pValue'};

% from the perspective of the depth 
tmp_depthR = [tbl.depth_right];
tmp_depthR(tmp_depthR == 0) = NaN; 
Mdl.depthR = fitlm([tbl.age], sqrt(tmp_depthR), 'VarNames', {'Age' , 'DepthR'}, 'RobustOpts', 'on');
handles4g = plotAdjustedResponse(Mdl.depthR, 'Age'); title ('Depth to younger vs. age');
adjusted4g = prepare_handle_for_python(handles4g);
adjusted4g.pval = Mdl.depthR.Coefficients{'Age', 'pValue'};

% from the perspective of the width 
tmp_widthR = [tbl.width_right];
tmp_widthR(tmp_widthR == 0) = NaN; 
Mdl.widthR = fitlm([tbl.age], tmp_widthR, 'VarNames', {'Age' , 'WidthR'}, 'RobustOpts', 'on');
handles4i = plotAdjustedResponse(Mdl.widthL, 'Age'); title ('Slope to younger vs. age');
adjusted4i = prepare_handle_for_python(handles4i);
adjusted4i.pval = Mdl.widthR.Coefficients{'Age', 'pValue'};


%% Where is the unstable point? Left or right? [not plotted]

% idx_left_tmp = tbl.depth_left>0;
% idx_right_tmp = tbl.depth_right>0;
% idx_both = [idx_right_tmp&idx_left_tmp];
% 
% idx_left = [idx_left_tmp&~idx_both]; 
% idx_right = [idx_right_tmp&~idx_both];
% 
% 
% [H3,P3,~,STATS3] = ttest2(tbl.age(idx_left), tbl.age(idx_right));

clearvars -except Mdl STATS2 STATS3 P3 P2 Val2_0 Val2_1 tbl

save('for_figure_3_landscapes_v9', 'n_unstable_points', 'age')
save('for_figure_4_landscapes_v9')


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
