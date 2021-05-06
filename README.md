# LinkHD


a versatile approach to integrate heterogeneous datasets



## Laura Zingaretti & Yuliaxis Ramayo Caldas

#### m.lau.zingaretti@gmail.com

If you find our package useful, please cite:

Laura M Zingaretti, Gilles Renand, Diego P Morgavi, Yuliaxis Ramayo-Caldas, Link-HD: a versatile framework to explore and integrate heterogeneous microbial communities, Bioinformatics, https://doi.org/10.1093/bioinformatics/btz862

LinkHD is a general R software to integrate heterogeneous dataset focusing on micribial communities. LinkHD combines multivariate techniques to perform data integration with cluster and variable selection.
The method also allows us to study the relashionships between observations and features and to obtain **enrichment taxa analysis**.

## Installation:

Clone the repository or
```{r}
devtools::install_github(repo="lauzingaretti/LinkHD")
```
or from Bioconductor

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 #The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("LinkHD")
```



## Usage
```{r}
library(linkHD)
```

LinkHD capabilities were demostrated analizing public datasets from TARA Ocean expedition (https://oceans.taraexpeditions.org/en/m/about-tara/les-expeditions/tara-oceans/) and datasets from rumen metataxonomic communities (including bactera, archaea and protozoa data). Data can be loaded by the next command:

```{r}
data(Rumynotipes)
#or
data(Tataoceans)
```
More examples and explanation of methods are available at:  https://lauzingaretti.github.io/LinkHD-examples/  


## Missing data

we have added a function to impute missing values in raw data. Note that this function should be used before any data transformation!

In fact, the process should be the following:


```{r}
#dir where data are stored
setwd('dir')
Data<-ReadData()
Out<-Na_inspect(Data)
# if some of your data contains NA, you can use impute_missing() for each data.frame, otherwise you can follow with the standard analysis.
```
