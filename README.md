![Build Status](https://travis-ci.com/JEstabrook/regulon-enrichment.svg?token=ZRDWBWe9sXCivP1NrZwq&branch=master)  [![](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-367) ![t](https://img.shields.io/badge/license-MIT-nrightgreen.svg) ![t](https://img.shields.io/badge/status-stable-nrightgreen.svg) ![t](https://zenodo.org/badge/179752059.svg)

# Priori

Priori is a Python module used to predict the activity of regulatory proteins from RNAseq data.

Priori submodules:

### `enricher.features` ###
Load -omic datasets


### `enricher.regulon` ###
Regulon utilities

# Dependencies

**Priori** requires:
~~~~~~~~~~~~
- Python (>= 3.6)
- scikit-learn (>= 0.21.3)
- NumPy (>= 1.17.3)
- SciPy (>= 1.3.1)
- pandas (>= 0.25.3)
- tqdm (>= 4.38.0)
- dill (>= 0.3.1.1)
~~~~~~~~~~~~

# Overview

Priori leverages pathway information and gene expression data to produce regulon-based protein activity scores. 
Our method tests for positional shifts in experimental-evidence supported networks consisting of transcription factors 
and their downstream signaling pathways when projected onto a rank-sorted gene-expression signature. 

# Running Priori

## Invoking Priori from the command line

Initialize github in the directory where you want to download Priori. Clone the Priori Github folder using
```
git clone https://github.com/ohsu-comp-bio/regulon-enrichment.git
```

Open the **regulon_enrichemnt** folder. Create a conda environment with the dependencies needed to run Priori
```
conda create -f priori_env.yml
```

Once the environment has been built, activate it
```
conda activate priori_env
```

Open the **enricher** folder. Set this path to your PATH variable. After sourcing your bashrc script, you should be able to run Priori using ::
    enrich

# Priori parameters

## Required parameters

`expr` : which tab delimited expression matrix to use shape : `[n_features, n_samples]`, units : `TPM, RPKM`

`out_dir` : output directory - directory serialized Enrichment object and enrichment.tsv will be saved to


## Optional parameters

`regulon` : optional regulon containing weight interactions between regulator and 
            downstream members of its regulon shape : `[len(Target), ['Regulator','Target','MoA','likelihood']`

`regulon_size` : number of downstream interactions required for a given regulator in order to calculate enrichment score `default=15`

`sec_intx` : path to pre-compiled serialized secondary interaction network, `default=secondary_intx_regulon.pkl`

`scaler_type` : scaler to normalized features/samples by: `standard | robust | minmax | quant`, default=`robust`

`thresh_filter` : Prior to normalization remove features that have a standard deviation per feature less than `{thresh_filter}`, `default=0.1`)


# Computing regulon enrichment scores

To quantify the regulon enrichment for a given dataset, the command line script `enrich` is used.

Use --help argument to view options

`enrich --help`

Priori requires two positional arguments: `expr` and `out_dir`

`enrich expr out_dir [regulon] [regulon_size] [sec_intx] [scaler_type] [thresh_filter] ` 

It is recommended to run enrich with the default parameters. 

`enrich tests/resources/test_expr.tsv test_enrichment_scores`

The command above will generate enrichment scores for the unittest dataset `test_expr.tsv` and will generate and store the output under `test_enrichment_scores/`. In this directory `test_enrichment_scores/`, both the serialized Enrichment object `test_enrichment.pkl` and a tsv of the enrichment scores,`test_regulon_enrichment.tsv` will be found. 

The `enrichment.tsv` file be shaped : `[n_samples, n_regulators]`, where `n_samples` refers to the original number of samples provided in `expr`, while `n_regulators` will be determined based on the overlapping features present in the `expr` dataset and the `regulon_size` parameter. 
