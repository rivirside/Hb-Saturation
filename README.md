# Hemoglobin Oxygen Saturation PMF Analysis

Interactive web-based simulator and analysis tools for exploring hemoglobin oxygen-binding probability mass functions (PMF) using binomial and Adair models.

## Overview

This repository contains tools for analyzing the occupancy distribution of hemoglobin tetramers, comparing independent-site (binomial) and cooperative (Adair) binding models. It accompanies the paper "Occupancy Distribution of Hemoglobin Tetramers: Analytic Probability-Mass Function and Physiological Implications" by Tomer Rivir Yonatan Zilbershtein.

## Features

### Interactive Web Simulator
- **Real-time visualization** of probability mass functions for both models
- **Adjustable pO₂** from 0-150 mmHg via slider control
- **Customizable K-values** (K₁-K₄) for the Adair model
- **Preset configurations** including Winslow constants
- **Side-by-side comparison** of binomial vs. Adair predictions

### Python Analysis Tools
- Comprehensive PMF calculations for both models
- Generation of publication-quality figures
- Sensitivity analysis of Adair binding constants
- RMSE validation against experimental data
- Export to CSV and LaTeX formats

## Quick Start

### Web Simulator
Simply open `index.html` in a modern web browser or visit the [live demo](https://rivirside.github.io/Hb-Saturation/).

### Python Analysis
```bash
# Install dependencies
pip install numpy pandas matplotlib scipy

# Run the analysis
python hb_pmf_analysis.py
```

## Models Implemented

### Binomial Model
- Assumes 4 independent oxygen-binding sites
- Uses Hill equation with P₅₀ = 26 mmHg, h = 2.7
- Probability: P(k) = C(4,k) × p^k × (1-p)^(4-k)

### Adair Model
- Accounts for cooperative binding
- Uses sequential binding constants K₁-K₄
- Default Winslow values: K₁=0.004, K₂=0.043, K₃=0.262, K₄=0.039 mmHg⁻¹

## Repository Structure
- `index.html` - Interactive web simulator
- `hb_pmf_analysis.py` - Python analysis script
- `LICENSE` - License information

## Citation
If you use these tools in your research, please cite:
```
Zilbershtein, T.R.Y. (2024). Occupancy Distribution of Hemoglobin Tetramers: 
Analytic Probability-Mass Function and Physiological Implications.
```

## License
This project is licensed under the terms in the LICENSE file.