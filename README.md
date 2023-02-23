# Priori

Priori predicts transcription factor activity from RNA sequencing data using prior, literature-supported regulatory relationship information.

# Running Priori

## Invoking Priori from the command line

Download the latest release using conda:
```
conda install -c wmyashar priori
```

Priori can also be downloaded directly from GitHub. 
```
git clone https://github.com/ohsu-comp-bio/regulon-enrichment.git
```

Open the **regulon-enrichment** folder. Create a conda environment with the dependencies needed to run Priori
```
cd regulon-enrichment
conda env create -f priori_env.yml
conda activate priori_env
```

Once the environment has been built, add the path to the **priori** folder to your PATH variable. After sourcing your bashrc script, you should be able to run Priori using
```
priori
```

# Priori parameters

## Required parameters

`expr` : which tab delimited expression matrix to use shape : `[n_features, n_samples]`, units : `TPM, RPKM`

`out_dir` : output directory - directory serialized Priori object and priori_activity_scores.tsv will be saved to


## Optional parameters

`regulon` : optional regulon containing weight interactions between regulator and 
            downstream members of its regulon shape : `[len(Target), ['Regulator','Target','MoA','likelihood']`

`regulon_size` : number of downstream interactions required for a given regulator in order to calculate enrichment score `default=15`

`sec_intx` : path to pre-compiled serialized secondary interaction network, `default=secondary_intx_regulon.pkl`

`scaler_type` : scaler to normalized features/samples by: `standard | robust | minmax | quant`, default=`robust`

`thresh_filter` : Prior to normalization remove features that have a standard deviation per feature less than `{thresh_filter}`, `default=0.1`)


# Computing transcription factor activity scores

To quantify the transcription factor activity for a given dataset, the command line script `priori` is used.

Use --help argument to view options

`priori --help`

Priori requires two positional arguments: `expr` and `out_dir`

`priori expr out_dir [regulon] [regulon_size] [sec_intx] [scaler_type] [thresh_filter] ` 
It is recommended to run Priori with the default parameters. 

The `priori_activity_scores.tsv` file be shaped : `[n_samples, n_regulators]`, where `n_samples` refers to the original number of samples provided in `expr`, while `n_regulators` will be determined based on the overlapping features present in the `expr` dataset and the `regulon_size` parameter. 
