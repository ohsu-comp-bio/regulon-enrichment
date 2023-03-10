# Priori

Priori predicts transcription factor activity from RNA sequencing data using prior, literature-supported regulatory relationship information.

# Configure

The latest release of Priori can be downloaded using conda:
```
conda install -c conda-forge priori
```

Alternatively, Priori can be downloaded directly from GitHub. 
```
git clone https://github.com/ohsu-comp-bio/regulon-enrichment.git
```

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

Estabrook, J, Yashar, WM, et al. Predicting transcription factor activity using prior biological information. bioRxiv (2023). https://doi.org/10.1101/2022.12.16.520295 
