# Introduction

Disruption of normal transcription factor regulation is associated with a broad range of diseases. It is important to detect aberrant transcription factor activity to better understand disease pathogenesis. We have developed Priori, a method to predict transcription factor activity from RNA sequencing data. Priori has two key advantages over existing methods. First, Priori utilizes literature-supported regulatory information to identify transcription factor-target relationships. It then applies linear models to determine the impact of transcription factor regulation on the expression of its target genes. Second, results from a third-party benchmarking pipeline reveals that Priori detects aberrant activity from 124 gene perturbation experiments with higher sensitivity and specificity than 11 other methods.

# Tutorial

We have created a tutorial that demonstrates how to generate and analyze Priori transcription factor activity scores. You can find the instructions to set up the environments and the Jupyter notebook demonstrating the use of Priori scores in the "tutorial" folder.

# Installation

The latest release of Priori can be downloaded using conda:
```
conda install -c conda-forge priori
```

Using pip:
```
pip install priori_regulon-enrichment
```

Or directly from GitHub: 
```
git clone https://github.com/ohsu-comp-bio/regulon-enrichment.git
```

# Set-up

Once Priori is downloaded, the input data needs to be formatted for analysis. Activity scores should only be generated from normalized gene expression data. Any standard normalization method, including CPM or TPM, stored in a tab-delimited file is sufficient. The input data must be in the following format:
  1. Gene symbols stored in a column labeled "features". 
  2. Separate columns of normalized gene expression values for each sample. Label each column with the sample name.

| features  | sample_1 | sample_2 | ... | sample_n |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| gene_1  | #  | # | ... | # |
| gene_2  | #  | # | ... | # |
| ...  | ...  | ... | ... | ... |
| gene_n  | #  | # | ... | # |

# Usage
```

priori expr out_dir [--help] [--regulon "<value>"] [--regulon_size "<value>"] 
                    [--scaler_type "<value>"] [--thresh_filter "<value>"] 

Required arguments:
    expr                A tab-delimited normalized expression matrix of the shape 
                        [n_features, n_samples]
                        
    out_dir             Output directory where the serialized Priori object and 
                        priori activity scores will be saved

Optional arguments:

    --regulon           A prior network that contains the transcriptional regulators 
                        (regulator), target genes (target), edge weights (moa), and
                        likelihood of interaction (likelihood). The network should be 
                        formatted as ['regulator','target','moa','likelihood']
                        
    --regulon_size      Number of downstream target genes required for a given 
                        transcriptional regulator. Default = 15
                        
    --scaler_type       Method to scale normalized expression. Options include standard, 
                        robust, minmax, or quant. Default = robust.
                        
    --thresh_filter     Remove features with a standard deviation below this value. Default = 0.1.
```

# Output

Priori generates three files: 
  1. `priori_activity_scores.tsv`: This file contains the normalized, transcription factor activity scores in the structure of [n_samples, n_features]. This is the primary output to analyze. Please see our [paper](https://doi.org/10.1016/j.isci.2024.109124) or the [tutorial](https://github.com/ohsu-comp-bio/regulon-enrichment/tree/master/tutorial) for ways to analyze these scores.
  2. `priori_activity_score_weights.tsv`: This file contains the transcription-factor target gene weights. This file is useful to understand the impact of individual target genes on the overall transcription factor activity scores. Target genes with a large absolute `MoA` value (relative to other target genes) have a greater impact on the activity score. A positive `MoA` indicates that the expression of a given target gene is directly correlated to the expression of the transcription factor (and vice versa).
  3. `priori_object.pkl`: This pickle stores all of the information that Priori used to generate the activity scores, including the input parameters. In order to load this file, use the `read_pickle` command from `regulon/regulon_utils.py`. 

# Paper

Priori has been published in Cell iScience. If you use our program in your studies, please cite our paper:

Yashar WM, Estabrook J, Holly HD, Somers J, Nikolova O, Babur Ö, Braun TP, Demir E. Predicting transcription factor activity using prior biological information. iScience. 2024 Feb 5;27(3):109124. doi: 10.1016/j.isci.2024.109124. PMID: 38455978; PMCID: PMC10918219.
