Use of control-augmentation to alleviate issues with cross-study transfer

Structure:
```
    -- files             results of cross predictions and such
    -- models            SIAMCAT models trained with control augmentation
    -- src               folder with the scripts needed to perform the analysis
    -- submission_temp   folder holding submission scripts
```

There are several scripts in this folder:

- `study_similarity.R`: This script tries to learn a simple ML model to
distinguish between control samples from different datasets (this will be used
    as basis for later scripts).
- `prep_ctr_augm_jobs.R` and `train_control_augmented.R`: These scripts
prepare and execute jobs on the cluster to train control-augmented models
(using the `utils.R` script), which come in different flavours:
    - `cohort_2`: control samples are added from cohort studies, n=2
    - `cohort_5`: control samples are added from cohort studies, n=5
    - `other_ctr`: control samples are added from other bigger studies
    (included in the ML meta-analysis)
    - `random`: whole studies from the ML meta-analysis are randomly added
    (including cases)
    - `similar`: studies are trained with control samples from similar
    studies added (based on `study_similarity`, see `control_augm_similar.R`)
- `cross_prediction.R`: how does control-augmentation change the
cross-prediction results for the models? The results are plotted by the script
`plot_cross_predictions.R`
- `weight_pcoa.R`: compare and plot the model weights before and after
control-augmentation
