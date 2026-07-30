# ExposureGAT: An Exposure-Aware Graph Attention Framework for Integrative Multi-Omics Data Analysis 


## Abstract
Multi-omics profiling can measure disease-related changes across gene expression, DNA methylation, copy-number variation, and other molecular data types. However, many existing approaches analyze these data as separate hypotheses or treat them as a flat feature matrix. They also often model exposure as a simple covariate or group label, limiting their ability to study how external exposures shape molecular changes between paired biological states. We developed ExposureGAT, an exposure-aware graph attention framework for paired multi-omics data analysis. ExposureGAT represents each participant using two matched graphs: one for a reference state and one for a transition state. In each graph, nodes represent molecular features, such as gene expression, DNA methylation, and copy-number variation, while edges represent known biological relationships, including gene-gene interactions and regulatory links. Unlike standard case-control classifiers, ExposureGAT learns a patient-level transition score that reflects how strongly paired molecular changes are associated with exposure intensity. The model also generates cohort-level and patient-specific node and edge scores, allowing important molecular features and rewired biological relationships to be identified across the cohort and within individual patients. We evaluated ExposureGAT using simulated paired multi-omics data and TCGA paired lung cancer data. Overall, ExposureGAT provides a flexible framework for studying how exposures are associated with shared and patient-specific molecular transitions across paired biological states.

## Usage 

```
ExposureGAT.py

Training script for ExposureGAT, an exposure-aware graph attention framework
for paired multi-omics graph learning.

This script trains ExposureGAT using paired reference-transition graphs and an
exposure-ranking objective. The model learns patient-level transition scores
associated with exposure intensity and saves the best model based on validation
Spearman correlation.

Example:
    python src/ExposureGAT.py \
      --data_dir GNN_Input_Data \
      --out_dir GNN_Output \
      --epochs 100
```

```
Analyze.py

Post-training analysis script for ExposureGAT.

This script loads a trained ExposureGAT model and extracts patient-level
transition scores, cohort-level node and edge scores, patient-specific node
scores, and rewiring outputs for downstream biological interpretation.

Example:
    python src/Analyze.py \
      --data_dir GNN_Input_Data \
      --model_path GNN_Output/best_model.pth \
      --config GNN_Output/run_config.json \
      --out_dir GNN_Output/results

```
### Key command-line arguments

| Argument | Description |
|---|---|
| `--graph_dir` | Directory containing paired graph input files |
| `--metadata` | Sample or patient metadata file |
| `--exposure_col` | Column containing the exposure variable |
| `--out_dir` | Output directory |
| `--epochs` | Number of training epochs |
| `--batch_size` | Training batch size |
| `--lr` | Learning rate |
| `--weight_decay` | Weight decay for optimizer |
| `--hidden_dim` | Hidden dimension for graph attention layers |
| `--dropout` | Dropout rate |
| `--seed` | Random seed |

