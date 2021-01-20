This folder is expected to contain all the cleaned tables from the
meta-analysis of
[Duvallet et al. _Nat Commun_ 2017](https://www.ncbi.nlm.nih.gov/pubmed/29209090).

Please go to
[this Github repository](https://github.com/cduvallet/microbiomeHD),
which collects all analysis scripts and results for the meta-analysis of
Duvallet et al.

After executing the `MAKEFILE` of the repository, the repository should contain
the cleaned abundance tables and metadata for all included datasets, which are
saved under
```
./data/clean_tables
ls | head
art_scher.metadata.clean.feather
art_scher.otu_table.clean.feather
asd_kang.metadata.clean.feather
asd_kang.otu_table.clean.feather
asd_son.metadata.clean.feather
asd_son.otu_table.clean.feather
cdi_schubert.metadata.clean.feather
cdi_schubert.otu_table.clean.feather
cdi_vincent.metadata.clean.feather
cdi_vincent.otu_table.clean.feather
```

Copy all the `.feather` files into this folder to perform
* the reproduction of the results from Duvallet et al. using SIAMCAT
* the parameter space exploration including the Duvallet datasets
