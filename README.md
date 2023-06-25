# cfDNA_classification_analyses

Repository containing the machine learning analyses based on the nucleosome positioning of cell-free DNA samples around binding sites of regulatory elements. These analyses were performed as part of a master thesis.

## About the features

Example feature data sets from own analyses is given in this folder. These features were extracted from scores calculated in the https://github.com/leaburkard/cfDNA_midpoint_coverage repository. As samples the data of Christiano *et al.* (Nature, 2019) was used and as regions the transcripion factor binding sites of the GTRD (version 19.10).

The features are compared to the feature set of the Griffin analyses published by Doebley *et al.* (Nature Communications, 2022). Their data is provided at https://github.com/adoebley/Griffin_analyses/tree/main/delfi_data_cancer_detection/number_of_sites_analysis/merged_data and can be downloaded from there.


## Generate own features
Execute the cfDNA_midpont_coverage repository with own cfDNA samples and regions of interest. Then use the feature_extraction.py script. 

```
-r: run ID given in the config/samples.tsv file in the cfDNA_midpont_coverage repository
-s: which of the three scores (MIDPOINT, WPS or COV) was calculated
-a: which method to calculate the amplitude (FFT or spec_pgram)
-p: path to the cfDNA_midpont_coverage repository
```

Example:
```
scripts/feature_extraction -r <ID> -s MIDPOINT -a FFT -p ../cfDNA_midpont_coverage
```

More information is given using the -h help option.

## Notebooks

The notebooks contain different machine learning analyses conducted with the described features.