# scRNA_Cluster

## What does this project do?
Recent technological advances now allow us to profile the transcriptome of individual cells, giving us insight into the RNA molecules present on a cell-by-cell basis. The amount of RNA corresponding to a given gene in a cell is indicative of the level of expression of that gene. In contrast to bulk RNA sequencing, which averages gene expression levels over thousands or even millions of cells, the finer resolution offered by single cell RNA (scRNA) data opens up many new avenues for exploration. 

The codes in this repository use data science and machine learning techniques to create a pipeline to cluster scRNA data from across 33 datasets from published papers. 

## Why is this project useful?
scRNA has been used to characterize new cell types and new cell states, including new lung cells and new brain cells. The more granular information also enables us to further our understanding of cell state transitions. scRNA profiles of individual cells can also be used to retrieve similar cell types. As an example application, suppose a population of cells is taken from a tumor. One question of interest is which immune cells are present in the tissue sample. Such analysis is often performed using known genetic markers, but a more comprehensive solution would compare the scRNA profile with a set of curated single cells with known types.

Below is a picture of an elbow plot by scripts in this code:
<img src="https://user-images.githubusercontent.com/26308648/48154913-be925400-e297-11e8-9c3c-a3686ddc5d36.png" width="620">

Example of one silhoutte plot created by this code can be found below:
<img src="https://user-images.githubusercontent.com/26308648/48155065-2ba5e980-e298-11e8-9359-68a157428ac7.png" width="620">

## How to get started with this project?
```
$ git clone https://github.com/sidd0529/scRNA_Cluster.git
$ cd scRNA_Cluster
```

Run the code using:
```
$ python cellpred.py
```


## Where can you get help with this project?
I will be very happy to help in case you have any questions regarding this project. You can find me at siddharthsatpathy.ss@gmail.com .
