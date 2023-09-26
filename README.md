# Priori

Disruption of normal transcription factor regulation is associated with a broad range of diseases. It is important to detect aberrant transcription factor activity to better understand disease pathogenesis. We have developed Priori, a method to predict transcription factor activity from RNA sequencing data. Priori has two key advantages over existing methods. First, Priori utilizes literature-supported regulatory information to identify transcription factor-target relationships. It then applies linear models to determine the impact of transcription factor regulation on the expression of its target genes. Second, results from a third-party benchmarking pipeline reveals that Priori detects aberrant activity from 124 gene perturbation experiments with higher sensitivity and specificity than 11 other methods.

# Tutorial

We have created a tutorial that demonstrates how to generate and analyze Priori transcription factor activity scores. You can find the instructions to set up the environments and the Jupyter notebook demonstrating the use of Priori scores in the "tutorial" folder.

# Installation

The latest release of Priori can be downloaded using conda:
```
conda install -c conda-forge priori
```

Alternatively, Priori can be downloaded directly from GitHub: 
```
git clone https://github.com/ohsu-comp-bio/regulon-enrichment.git
```

# Set-up

Once Priori is downloaded, the input data needs to be formatted for analysis. Activity scores should only be generated from normalized gene expression data. Any standard normalization method, including CPM or TPM, stored in a tab-delimited file is sufficient. The input data must be in the following format:
  1. Gene symbols stored in a column labeled "features".
  2. Separate columns of normalized gene expression values for each sample. Label each column with the sample name.

| features  | Sample_1 | Sample_1 | ... | Sample_n |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Gene_1  | #  | # | ... | # |
| Gene_2  | #  | # | ... | # |
| ...  | #  | # | ... | ... |
| Gene_n  | #  | # | ... | # |

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
                        (Regulator), target genes (Target), edge weights (MoA), and
                        likelihood of interaction (likelihood). The network should be 
                        formatted as ['Regulator','Target','MoA','likelihood']
                        
    --regulon_size      Number of downstream target genes required for a given 
                        transcriptional regulator. Default = 15
                        
    --scaler_type       Method to scale normalized expression. Options include standard, 
                        robust, minmax, or quant. Default = robust.
                        
    --thresh_filter     Remove features with a standard deviation below this value. Default = 0.1.
```

# Paper

Priori has been released as a pre-print. If you use our program in your studies, please cite our paper:

Yashar, WM, Estabrook, J, et al. Predicting transcription factor activity using prior biological information. bioRxiv (2023). https://doi.org/10.1101/2022.12.16.520295 
