![Build Status](https://travis-ci.com/JEstabrook/regulon-enrichment.svg?token=ZRDWBWe9sXCivP1NrZwq&branch=master)  [![](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-367) ![t](https://img.shields.io/badge/license-MIT-nrightgreen.svg) ![t](https://img.shields.io/badge/status-stable-nrightgreen.svg)


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
the easiest way to install regulon-enrichment is using ``pip``   ::

    pip install regulon-enrichment==0.0.2b0

or ``conda``::

    conda install -c estabroj89 regulon-enrichment


# Overview

This method leverages pathway information and gene expression data to produce regulon-based protein activity scores. 
Our method tests for positional shifts in experimental-evidence supported networks consisting of transcription factors 
and their downstream signaling pathways when projected onto a rank-sorted gene-expression signature. A regulatory gene’s
expression is generally a poor proxy for its protein’s behavior. We shift our focus to the collective enrichment of a 
regulatory protein’s target genes (its “regulon”). We extract protein-protein edges from Pathway Commons to build 
regulatory networks. 

