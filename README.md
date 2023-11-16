# Identifying important pairwise logratios (PLRs) using sparse PCA (sPCA)
The proposed method to identify important pairwise logratios (and parts) using sparse PCA (introduced in "Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N., Aravkin, A. Y. (2020). Sparse principal component analysis via variable projection. SIAM Journal on Applied Mathematics 80(2):977–1002") is introduced in "Nesrstová, V., Wilms, I., Hron, K., Filzmoser, P. Identifying Important Pairwise Logratios in Compositional Data with Sparse Principal Component Analysis" (currently under revision; November 2023).

## Requirements
Calculations were performed in following software, using the listed libraries:
- R version: 4.1.3 (2022-03-10) -- "One Push-Up"
- Essential libraries:
    - sparsepca: for sparse PCA; version 0.1.2
    - easyCODA: for STEP; version 0.34.3
    - compositions: to handle compositional data; version 2.0-4

## Available scripts
1. Auxiliary functions
2. Simulations

   a) Scenario A

   b) Scenario B

   c) Scenario C
4. Real data applications

   a) Aar massif data set

   b) Archaeometric data set
6. Toy example

An easy-to-use example of how to perform the calculations and obtain the matrix of sparse loadings.
