# HumanDLPFC24
A cellular map of the DLPFC brain region and associations to Alzheimer's disease and cellular communities, from single nucleus RNA-seq of 
24 individuals and bulk RNA profiles of 640 individuals with estimated cellular compositions by the CelMod algorithm.
Cite: https://doi.org/10.1101/2020.12.22.424084

1. To run the mediation analysis linking cell populations to disease cascde use the following code and follow the documentation within the code:
scRNA_n24_YangHS_16AUG2021_edit29MAR2023.R

2. The code for the CelMod algorithm can be found here : https://github.com/MenonLab/Celmod

3. To run the cellular community analysis run the following:
Use the example input that can be found in the compressed folder: input_for_frequencies.zip

There are 2 main functions: 
plot_frequencies_and_cor(dfs, group.by, ident.by, path)
For each cell type, the function plots the proportions of the subcluster (ident.by column) per individual (or any other category - group.by column) 
and the correlations of these frequencies. It plots the frequencies for each cell type separately and also for all the cell types together.
dfs is assumed to be a list of dataframes, such that each data frame has (at least) 3 columns: "cell.name"   group.by  ident.by.
See examples in the: input_for_frequencies 
For example, using a Seurat object, with a column in meta.data "individual": 
data.frame(cell.name = rownames(obj@meta.data), group.is =  obj@meta.data$individual, ident.by = Idents(obj))
That should create the dataframe to use as input for a single cell type

The second function plots the networks and uses the igraph package 
Regarding the arguments, it uses the same list of dataframes as the previous function but in addition needs a dataframe of information to use per cell subtype.
This  infos dataframe is used to "color" the nodes of the network
The rownames of the infos dataframe should match the names of the subtypes cluster uses in the list of dataframes - there is also an example in the input folder

The function has few other parameters: bg.by, frame.by, gray.by
The bg.by is the background color of the nodes. It has to be a column name from the infos dataframe. frame.by is the border color - can be empty.
The gray.by is if you want some nodes to be gray as "disable" - it's a vector of 2 values: the first is the column of infos to use. 
And the second is the value to put as gray (see the example at line 197 in the code (frequencies.R)).

And, lastly, if you want to calculate a community detection algorithm, use plot.community= T (for example line 190). 
For now it uses a leading eigenvector method but igraph has a lot of options you can use (change calculate.communities line 131 of network_construction.R). 

Also, you can control the threshold to plot the correlations or anti-correlations using thres.edge and thres.edge.neg - 
changing this threshold changes the structure of the network you get, so it could be important. 
The structure is built only using the positive correlations pairs. You also can input your own structure (as a dataframe of coordinates).
