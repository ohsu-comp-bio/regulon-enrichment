# Introduction

This tutorial will illustrate how to generate and analyze Priori transcription factor activity scores. For this demonstration, we will replicate the first steps of the TCGA BRCA analysis from our [pre-print](https://www.biorxiv.org/content/10.1101/2022.12.16.520295v2).

This notebook requires the use of the conda package manager. If you do not have anaconda or miniconda3 already installed on your computer, please first follow the steps outlined in [this helpful blogpost](https://engineeringfordatascience.com/posts/install_miniconda_from_the_command_line/) or from the [miniconda documentation](https://docs.conda.io/projects/miniconda/en/latest/) to set it up on your device.

# Set-up

Once conda has been set-up, we need to install the appropriate packages in order to run our analysis. First, download the files in this tutorial folder to your personal device, unzip it, and then nativigate into the folder:
```
wget --no-check-certificate https://github.com/ohsu-comp-bio/regulon-enrichment/raw/master/tutorial/tutorial.zip -O tutorial.zip
unzip tutorial.zip
cd tutorial
```

Navigate to the tutorial file in the command line and then run the command below:
```
conda env create -f priori_tutorial_env.yml
```

When this operation is finished, activate the environment with the following command:
```
conda activate priori_tutorial_env
```

You can now initialize jupyter notebooks from the command line and proceed with the tutorial in the document titled `priori_tcga_brca_tutorial.ipynb`:
```
jupyter notebook
```