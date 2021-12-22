![Build Status](https://travis-ci.com/JEstabrook/regulon-enrichment.svg?token=ZRDWBWe9sXCivP1NrZwq&branch=master)  [![](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-367) ![t](https://img.shields.io/badge/license-MIT-nrightgreen.svg) ![t](https://img.shields.io/badge/status-stable-nrightgreen.svg) ![t](https://zenodo.org/badge/179752059.svg)


# Enrich

**regulon-enrichment** is a Python module used to predict the activity of regulatory proteins from RNAseq data.

*regulon-enrichment* submodules:

### `enricher.features` ###
Load -omic datasets


### `enricher.regulon` ###
Regulon utilities

# Dependencies

**regulon-enrichment** requires:
~~~~~~~~~~~~
- Python (>= 3.6)
- scikit-learn (>= 0.21.3)
- NumPy (>= 1.17.3)
- SciPy (>= 1.3.1)
- pandas (>= 0.25.3)
- tqdm (>= 4.38.0)
- dill (>= 0.3.1.1)
~~~~~~~~~~~~

# User installation
~~~~~~~~~~~~~~~~~

If you already have a working installation of numpy and scipy,
the easiest way to install regulon-enrichment is using ``conda``   ::

    conda install -c estabroj89 regulon-enrichment

or ``pip``::

    pip install regulon-enrichment==0.0.2b0

~~~~~~~~~~~~~~~~~
# Overview

This method leverages pathway information and gene expression data to produce regulon-based protein activity scores. 
Our method tests for positional shifts in experimental-evidence supported networks consisting of transcription factors 
and their downstream signaling pathways when projected onto a rank-sorted gene-expression signature. 

This regulon enrichment method utilizes pathway and molecular interactions and mechanisms available through Pathway 
Commons to accurately infer aberrant transcription factor activity from gene expression data.

# Running regulon-enrichment
## Invoking enrich from the command line

When installing the regulon-enrichment package, the set of scripts that make up to inteface to regulon-enrichment will 
automatically be placed as an executables in your path, so that you can refer to these without modifying your shell 
environment. For example, if you install regulon-enrichment using conda, then enrich will become available on the path, 
and  can be run as:

~~~~~~~~~~~~~~~~~
enrich
~~~~~~~~~~~~~~~~~

# Enrich parameters

## Required parameters

`cohort` : which cohort to use; this information will be retained in the serialized Enrichment class

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

Enrich requires three positional arguments: `cohort`,`expr`, `out_dir`

`enrich cohort expr out_dir [regulon] [regulon_size] [sec_intx] [scaler_type] [thresh_filter] ` 

It is recommended to run enrich with the default parameters. 


`enrich test tests/resources/test_expr.tsv test_enrichment_scores`

The command above will generate enrichment scores for the unittest dataset `test_expr.tsv` and will generate and store the output under `test_enrichment_scores/`. In this directory `test_enrichment_scores/`, both the serialized Enrichment object `test_enrichment.pkl` and a tsv of the enrichment scores,`test_regulon_enrichment.tsv` will be found. 

The `enrichment.tsv` file be shaped : `[n_samples, n_regulators]`, where `n_samples` refers to the original number of samples provided in `expr`, while `n_regulators` will be determined based on the overlapping features present in the `expr` dataset and the `regulon_size` parameter. 
