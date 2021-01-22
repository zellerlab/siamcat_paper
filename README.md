# Code and Analysis Results for the SIAMCAT paper

This repository contains a collection of smaller (and larger) analyses that
are included in the SIAMCAT paper.

For more information, please see:

> Jakob Wirbel, Konrad Zych, Morgan Essex, Nicolai Karcher, Ece Kartal,
> Guillem Salazar, Peer Bork, Shinichi Sunagawa, Georg Zeller  
> **Microbiome meta-analysis and cross-disease comparison enabled by the
> SIAMCAT machine-learning toolbox**
> Preprint on https://www.biorxiv.org/content/10.1101/2020.02.06.931808v2 (2020)


##### Requirements

All analyses presented here were conducted in **R version 3.6.1** and with the
[SIAMCAT v1.5.1](https://github.com/zellerlab/siamcat/releases/tag/v1.5.1).

In order to reproduce the results, several R packages are needed. You can
check if you have all packages installed by running:
```bash
Rscript ./utils/test_requirements.R
```

The `ggembl` package is only used for plotting and can be found
[on the EMBL Gitlab system](https://git.embl.de/grp-zeller/ggembl).

##### Reproducibility

Please note that re-running the code might lead to results that are numerically
not exactly the same as the results reported in the paper. This is due to
random seeds affecting cross-validation splits.  
Qualitatively, though, the results should mirror what is reported in the paper.

##### Preparations

In order to prepare the analyses, several datasets have to be downloaded from
Zenodo and other places. You can start the process by running:
```bash
Rscript ./utils/prepare_data.R
```

##### Taxonomic and functional profiling

For the presented analysis, several metagenomic datasets have been processed
with [_ngless_](https://www.ncbi.nlm.nih.gov/pubmed/31159881) in order to
create taxonomic profiles with the
[mOTUs2](https://www.ncbi.nlm.nih.gov/pubmed/30833550) tool and functional
profiles by mapping read to
[eggNOG4.5](https://www.ncbi.nlm.nih.gov/pubmed/26582926) groups.

The script to process the datasets can be found under
[ngless_script](https://github.com/zellerlab/siamcat_paper/tree/master/utils/full_pipeline.ngl).


## Primary outputs

>See [primary_outputs](https://github.com/zellerlab/siamcat_paper/tree/master/primary_outputs)

To illustrate the primary outputs of the package, we used SIAMCAT to
analysed the data from
[Nielsen et al. _Nat Biotechnol_ 2014](https://www.ncbi.nlm.nih.gov/pubmed/24997787)
which are available through the
[_curatedMetagenomcisData_ R package](https://waldronlab.io/curatedMetagenomicData/).

This analysis corresponds to **Figure 1** in the manuscript.

## Identification of confounders

>See [primary_outputs](https://github.com/zellerlab/siamcat_paper/tree/master/primary_outputs)
and [metformin](https://github.com/zellerlab/siamcat_paper/tree/master/metformin)

To show how SIAMCAT can aid to detect confounders in metagenomic studies,
we use two examples:
1. The data from
[Nielsen et al. _Nat Biotechnol_ 2014](https://www.ncbi.nlm.nih.gov/pubmed/24997787)
contain control samples from both Spain and Denmark, but disease samples were
taken only from Spanish individuals. We show how the factor `Country` could
bias the results of a naive analysis.
2. Using the data from
[Forslund et al. _Nature_ 2015](https://www.ncbi.nlm.nih.gov/pubmed/26633628),
we used SIAMCAT to reproduce the analyses that showed a big influence of
metformin treatment on the microbiome composition of Type 2 diabetes cases.

These analyses corresponds to **Figure 2** and **Supplementary Figure 2**
in the manuscript.


## Illustration of machine learning pitfalls

>See [ml_pitfalls](https://github.com/zellerlab/siamcat_paper/tree/master/ml_pitfalls)

In this analysis, we demonstrate how two common machine learning pitfalls,
i.e. supervised feature selection and naive splitting of dependent data, can
lead to over-optimistic estimations of performance.

This analysis corresponds to **Figure 3** in the manuscript.

## Parameter space exploration

>See [parameter_space](https://github.com/zellerlab/siamcat_paper/tree/master/parameter_space)

In order to give meaningful recommendations about parameter choices with the
machine learning pipeline, we analysed a large set of case-control metagenomic
datasets that had been processed with a wide variety of taxonomic and
functional profiling tools.

This analysis corresponds to **Figure 4** and **Supplementary Figures 4-8**
in the manuscript.

## Control-augmentation

>See [parameter_space](https://github.com/zellerlab/siamcat_paper/tree/master/control_augmentation)

After we discovered dramatic problems when transferring naive ML models across
datasets, we devised a control-augmentation strategy, which works by adding
external controls during ML model training.

This analysis corresponds to **Figure 5** and **Supplementary Figures 9-14**
in the manuscript.

## IBD meta-analysis

>See [ibd_meta_analysis](https://github.com/zellerlab/siamcat_paper/tree/master/ibd_meta_analysis)

To demonstrate how SIAMCAT can enable machine learning meta-analyses of
metagenomic datasets, we performed a meta-anlaysis of five different
inflammatory bowel disease (IBD) metagenomic studies.

This analysis corresponds to **Figure 6** in the manuscript.

## Other analyses

**1. Ocean metagenomic analysis:** SIAMCAT can not only be used for human gut
metagenomics samples, but also for environmental and/or metatranscriptomic
samples. We show this by analysing data from
[Salazar et al. _Cell_ 2019](https://www.ncbi.nlm.nih.gov/pubmed/31730850).
This analysis corresponds to **Supplementary Figure 15** in the manuscript.
>See [ocean](https://github.com/zellerlab/siamcat_paper/tree/master/ocean)


**2. Reproduction of previous machine learning analyses:** We used SIAMCAT to
reproduce the results of two recent machine learning meta-analyses for
metagenomic data
([Pasolli et al. _PLoS Comput Biol_ 2016](https://www.ncbi.nlm.nih.gov/pubmed/27400279)
and
[Duvallet et al. _Nat Commun_ 2017](https://www.ncbi.nlm.nih.gov/pubmed/29209090)).
This analysis corresponds to **Supplementary Figure 1** in the manuscript.
>See [reproduction](https://github.com/zellerlab/siamcat_paper/tree/master/reproduction)

## Contact

For questions, suggestions, or in the case of any problems, please feel free
to contact me, [Jakob Wirbel](mailto:jakob.wirbel@embl.de).
