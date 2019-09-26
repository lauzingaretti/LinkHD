# LinkHD
a versatile approach to integrate heterogeneous datasets
## Laura Zingaretti & Yuliaxis Ramayo Caldas

#### m.lau.zingaretti@gmail.com

If you find this resource useful, please cite: LinkHD: a versatile framework to explore and integrate heterogeneous microbial comunnities. submitted. 

LinkHD is a general R software to integrate heterogeneous dataset focusing on micribial communities. LinkHD combines multivariate techniques to perform data integration with cluster and variable selection.
The method also allows us to study the relashionships between observations and features and to obtain **enrichment taxa analysis**. 

## Installation:

Clone the repository or 

devtools::install_github(repo="lauzingaretti/LinkHD")

## Usage

library(linkHD)

To illustrate the LinkHD capabilities we used a public data from TARA Ocean experidition as well as from the ruminal metataxonomic communities (including bactera, archaea and protozoa data). 
Both dataset are including in the package.

data(Rumynotipes)

data(Tataoceans)

Yo can found examples and methods explanation in the following web site: https://lauzingaretti.github.io/LinkHD/
