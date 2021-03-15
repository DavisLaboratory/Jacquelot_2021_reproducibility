This repository contains data and source code used to train the ILC2 XGBoost classifier in the manuscript ***PD-1 blockade unleashes ILC2-dependent anti-tumor immunity in melanoma*** by Jacquelot et al. The XGBoost model is provided as a `RData` object for easy application. Install `xgboost` package from CRAN in R >= 3.6.3, then:

### Load pre-trained XGBoost model

```
library(xgboost)
load("tonsils_ILC2_xgb_model.RData")

```

Access classifier feature (gene) names by
```
head(bst$feature_names)
nFeatures <- length(bst$feature_names)
```

### Apply the model to new datasets
To use the model to obtain ILC2 prediction scores on new datasets, apply TMM-normalisation to your RNAseq count data (optional), compute log2-cpms and subset your dataset on genes in the classifier (Note feature names are Homo-sapiens gene symbols). Here `testdata` is a matrix of logCPM values with `nFeatures` columns and `nSamples` rows

```
preds <- predict(bst, testdata)
```



