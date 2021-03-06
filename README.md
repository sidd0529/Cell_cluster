# Cell_cluster

## What does this project do?
Through recent technical advances, we can now profile the transcriptome of individual cells. Such research can help us look into RNA molecules inside cells. The amount of RNA in a given cell can give us knowledge of the level of expression of different genes in the cell. Such research offer us fin resolution into the expression of genes in millions of cells. 

The code in this repository uses data from 33 different labs and research publications and researches methods to predict cell types based on gene expression data. In addition to cell-gene data, this work also uses protein-DNA and protein-protein data to find similar genes. It also uses the said data to perform dimensionality reduction. The raw data has 17382 cells and over 80,000 features (genes). The pipeline that is presented in this repository uses novel methods for dimensionality reduction, good practices for data processing and tools in machine learning like Principal Component Analysis (PCA) and Logistic Regression to predict cell types. The prediction accuracy that is gotten in this code is close to 95%.

## Why is this project useful?
One can use the methods ilustrated in this code to identify new cell types (like new lung cells and new brain cells) and find which genes contribute to cell types.

Below is a picture of explained variance ratio obtained during implementation of PCA in this code:
<img src="https://user-images.githubusercontent.com/26308648/50857343-56292700-135b-11e9-80a8-491e7b48b691.png" width="620">

A plot which shows segregation of data (in the test data set) along the top two principal axes can be found below:
<img src="https://user-images.githubusercontent.com/26308648/50857361-66410680-135b-11e9-98fb-5152d7cfe632.png" width="620">

## How to get started with this project?
```
$ git clone https://github.com/sidd0529/Cell_cluster.git
$ cd Cell_cluster
```

One can download the folder containing files with cell-gene data from 33 labs [here](https://drive.google.com/open?id=1pBkCNAZ2eNSWFN5lcM7CaoBRqp9agg_E). Folder with protein-DNA interaction data can be gotten from [here](https://drive.google.com/open?id=16A92yjaU4vNsMTvk3_D3nyARjf3QI08c) . Download these folders and put them inside the folder 'Cell_cluster'.

Run the code using:
```
$ python cell_cluster.py
```


## Where can you get help with this project?
I will be very happy to help in case you have any questions regarding this project. You can find me at siddharthsatpathy.ss@gmail.com .
