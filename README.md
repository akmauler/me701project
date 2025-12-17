# me701project
# Combined Dosimetry Analysis and Statistical Error Propagation

This repository contains `ReadPHITSOutput.py`, a Python-based tool designed to automate the post-processing of Monte Carlo dosimetry data from PHITS (Particle and Heavy Ion Transport code System). The script is optimized for combining multiple radiation fields and performing rigorous statistical error propagation for radiotherapy research.

## Overview

In radiation dosimetry, accurate dose calculation is essential for understanding the biological and physical effects of absorbed radiation. While PHITS is the gold standard for modeling these interactions, manual summation of multiple beam outputs is tedious and error-prone. 

This script provides a standardized, traceable workflow for:
* **Dose Normalization**: Scaling raw dose output (Gy/source) to clinical dose (Gy) using field-specific Monitor Units (MU).
* **Automated Summation**: Combining independent radiation fields into a single total absorbed dose per anatomical region.
* **Error Propagation**: Utilizing the Root-Sum-Square (RSS) method to ensure accurate statistical confidence.

## Statistical Methodology

Because the simulations for each radiation field are statistically independent, a linear sum of errors would incorrectly overestimate the final uncertainty. The script calculates the total absolute error ($\sigma_{total,j}$) for any given region $j$ using the RSS method:

$$\sigma_{total,j}=\sqrt{\sum_{i=1}^{N}\sigma_{i,j}^{2}}$$

Where:
* $\sigma_{i,j}$ is the absolute error (standard deviation) of the scaled dose for field $i$ in region $j$.
* $N$ is the total number of simulated radiation fields.

## Key Features

* **Library Integration**: Leverages **Pandas** for efficient data handling, **NumPy** for vectorized RSS calculations, and **Matplotlib** for visualization.
* **Universal Parser**: Specifically designed to read PHITS `T-Deposit` tallies with `mesh = reg`, making it compatible with various human phantom models.
* **Error Handling**: Implements `try-except` blocks to prevent script failure due to missing files or improper formatting.

## Usage

### 1. Configuration
Update the `FILES_AND_MULTIPLIERS` dictionary in the script with your output file paths and their respective scaling factors:

```python
FILES_AND_MULTIPLIERS = {
    "outputfiles/Result_ATP-AM_VARIAN.out": 72,
    "outputfiles/Result_ATP-AM_VARIAN-field2.out": 71,
    # Add additional fields as needed
}
