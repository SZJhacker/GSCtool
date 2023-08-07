# GSCtool: a novel descriptor that characterizes genome for applying machine learning in genomics
Below is the sample data and code in the paper.
## 1. Sample data
**Sample data list**:
- dataA.gsc.csv: Feature extraction file for dataA.
- dataV.gsc.csv.gz: Feature extraction file for dataB, To use it, uncompress it by running `gunzip dataB.gsc.csc.gz`.
- dataC.gsc.csc.gz: Feature extraction file for dataC, To use it, uncompress it by running `gunzip dataC.gsc.csc.gz`.
- gl.csv: Grain length for all datasets
- gl.csv: Grain width for all datasets
  
**Description**
---
+ dataA.gsc.csv, dataB.gsc.csv and dataC.csv are GSCtool feature extracted from VCF files by GSCtool.
+ gw.csv and gl.csv are phenotype data files that provide phenotypic information related to GSCtool features.
  
## 2. Code
### 2.1 Run environment
#### Anaconda Installation
Make sure you have installed Anaconda or Miniconda, which will provide a unified Python environment manager and package manager. You can download and install Anaconda from the official Anaconda website (https://www.anaconda.com/download).

#### Dependency Package Installation
Install the following dependency packages to run the code:
```bash
conda install -c conda-forge pandas optuna scikit-learn lightgbm tensorflow
```

## 2.2 Code file
**DL_prediction.ipynb** - Jupyter Notebook file containing implementation and training code for Optuna and deep learning models.

Open the DL_prediction.ipynb using Jupyter Notebook and run each code block in the order it is listed. Make sure that the data files have been placed in the correct location and that the data file paths are specified correctly in the code.