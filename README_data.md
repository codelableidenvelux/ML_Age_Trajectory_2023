# Data Package for __Age-related behavioral resilience in smartphone touchscreen interaction dynamics__
Enea Ceolini, K. Richard Ridderinkhof, and Arko Ghosh

This data package contains 3 mat files:

Age trajectories predictions from BAM: all_age_pred_21d_30_11_2023.mat
    - all_ages: choronlogical age of each subject. 
    - all_genders: gender of each subject. 
    - all_pub: QuantActions' public ID.
    - all_preds: trajectories predicted from BAM.
    - all_preds_times: trajectories times (UTC).


Mean exit times for each subject/trajectory: all_mean_exit_times_28_12_2023.mat

Langevin analysis output: analysis_output_landscapes_15_12_2023.mat
    - all_mods: langevin models for each subject/trajectory
    - all_MU / all_SIGMA: for each subject/trajectory to reproject the analysis in the "age" space since the trajectories are first z-scored. 