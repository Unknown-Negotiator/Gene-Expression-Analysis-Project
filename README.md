# Gene-Expression-Analysis-Project
## [Project Code](https://colab.research.google.com/drive/1H5h8SE81VHLBebmaN5ZppS3MgYLbXnhS?usp=sharing)
## [Project Report](https://docs.google.com/document/d/1Z-6b-gMOwuYj3z23w6n1fa7wI4TbhQczHt4jumNr7ws/edit?usp=sharing)

# Exploratory Data Analysis
## Raw data distribution
Raw Counts |	Log2 Transformed Counts
-|-
![](images/raw_counts_histplot_2.png) | ![](images/log2_counts_histplot.png)

## Control vs Treated replicates boxplot
![](images/log2_counts_boxplot.png)
Here we can clearly observe a batch effect given by replicate: boxplots related to a similar replicate, but different in terms of treatment have more similar to each other expression distributions.

# Data Preprocessing
## Filtration
We filter out genes that have less than 10 read counts in total
Before filtration |	After filtration
-|-
![](images/counts_pre_10-filter_cut.png) | ![](images/counts_post_10-filter.png)

As a result of this filtration the **number of genes came down from 33602 to 24303**

## Normalization

# Differential Gene Expression Analysis
## MA plot
![](images/MA_plot.png)

## DEG plots
Before filtration |	After filtration
-|-
![](images/volcano_plot.png) | ![](images/control_vs_treatment.png)
