# Age-related behavioral resilience in smartphone touchscreen interaction dynamics    
Enea Ceolini, K. Richard Ridderinkhof, and Arko Ghosh

## Behavioral age model (BAM)
The training proceedure and output analysis of the BAM can be found in [this notebook](train_BAM.ipynb).

## Comparison bin method vs mesh method
We compare the output of Langevin reconstuction using two different methods: bin and mesh.
The analysis on some [sample data](examples_ages.mat) can be found [here](comparison_mesh_bins.m).

## Langevin reconstruction
We run the langevin reconstruction [here](full_langevin_analysis_mesh.m) on all the available age trajectories.

## Discern younger and older states
Younger and older states [can be identified](younger_older_bucket_states.m) in subjects with 2 stable points and 1 tipping point.

## Full analysis 
Fianl and full analysis of the results can be found [here](analysis_and_states.m).

## Bootstrap and null hypotesis testing
Following [this work](https://www.pnas.org/doi/10.1073/pnas.0802430105#supplementary-materials)
We want to see in what cases a surrgate time series presents alternative stable states by chance.
See [here](null_distributions.m).
 
for each sample:
    1. We create 1000 surrogate times series (each of 3 methods)
    2. We calculate the potential_eff
    3. We see what percentage of those have >= number of stable points w.r.t the original time series
    4. We see in how many of the samples the results has P. 0.05
Now for the 3 methods
    a. simple randomization with replacement 
    b. same spectrum
    c. autoregressive model  

## Acknowledgments
This work builds on top of the work done in:
- [Exit time as a measure of ecological resilience; Arani et. al](https://www.science.org/doi/10.1126/science.aay4895) 
- [Early Warnings of Regime Shift When the Ecosystem Structure Is Unknown
William A. Brock, Stephen R. Carpenter](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0045586)
- [A model of healthy aging based on smartphone
interactions reveals advanced behavioral age in
neurological disease; Ceolini et. al](https://doi.org/10.1016/j.isci.2022.104792)

## Miscellaneous 

Xgboost default parameters can be fiund on [their documentation](https://xgboost.readthedocs.io/en/release_1.7.0/parameter.html) 
The grid search was performed for the following parameters and corresponding values indicated in ‘[]’ – ‘max_depth’ [2, 4, 6, 8] ;  'min_child_weight' [2, 4, 6, 8] ; 'gamma' [0 to 0.5]; 'subsample' [0.6 to 1]; 'colsample_bytree' [0.1 to 0.7]; 'reg_alpha' [0, 0.001, 0.005, 0.01, 0.05]; and 'reg_lambda' [1x10-5, 1x10-2, 0.1, 1, 100].
