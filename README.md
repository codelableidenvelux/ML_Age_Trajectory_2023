# Age-related behavioral resilience in smartphone touchscreen interaction dynamics    

Enea Ceolini1, K. Richard Ridderinkhof, and Arko Ghosh


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
