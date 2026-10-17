# Environments

**QIIME 2** (steps 0–6 and 8) — official amplicon distribution, release 2024.10:

```bash
conda env create -n qiime2-amplicon-2024.10 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
```

**R** (step 7) — phyloseq, vegan, lme4/lmerTest:

```bash
conda env create -f envs/r-analysis.yml
```
