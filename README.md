# gene-expression-analysis

### Gene differential expression and disease prediction model

##### Data preparation:
Mice Protein Expression Data Set is downloaded from UCI’s Machine Learning Repository (https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression).

Here we focus on genes differentially expressed between control mice and trisomic mice. To avoid influence by the injection treatment, we omit the measurements belonging to the memantine category and retain the saline group. Moreover, measurements including NA value are excluded and finally we obtain a dataset with 297 measurements on 77 protein/protein modifications (150 for control mice and 147 for trisomic mice).

##### Method:
We implement **multiple t-test** to determine the statistical significance of differential expression for each protein/protein modification. To deal with the negative positive problem, we use the **Bonferroni correction** to adjust the p-value. Other methods applicable for adjusting p-value including **Benjamin-Hochberg (BH)** and **FDR**, are also available in our program. 

Based on the differentially expressed proteins/protein modifications, **LDA (Linear Discriminant Analysis)** and **10-fold cross-validation** are used to generate the predicion model.

##### Code:
R language is used and all involved source code is included in gene_diff.R.
