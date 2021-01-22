Exploration of the space of possible machine learning workflow hyper-parameters.

Structure:
```
    -- data_info         information about tasks
    -- files             miscellaneous useful information
    -- job_info          information about hyper-parameters
    -- models            SIAMCAT models trained with the optimal parameter set
    -- results           folder holding the result tables
    -- sc                folder containing the SIAMCAT objects
    -- src               folder with the scripts needed to perform the analysis
    -- submission_temp   folder holding submission scripts
    -- temp              folder to hold partial results
```

In short, the different scripts perform the following tasks:
1. `1_prep_datasets.R`: prepare the data and save `SIAMCAT` objects for each
task
2. `2_prep_parameter_space_jobs.R`: prepare parameter space for the exploration
3. `3_train_models.R`: train models on the EMBL cluster
4. `4_combine_results.R`: combine the partial results into a full matrix

The other scripts perform different analyses based on the results of the
parameter space exploration. They should be run more or less in this order
because some scripts might produce some files that later scripts need.

- `5_plot_results.R`: Plot the results of the large scale application
- `6_variance_explained.R`: which hyper-parameter in the machine learning
    workflow has the biggest influence?
- `7_permanova_auroc.R`: Compare the performance of common ecological
    distances and ML approaches to separate between cases and controls.
- `8_prep_naive_jobs.R` and `9_train_best_models.R`: Train models for the
    best-performing parameter set
- `10_cross_predictions.R`: How easily can ML models be transferred across
    studies?
