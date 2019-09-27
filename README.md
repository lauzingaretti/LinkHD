# LinkHD


a versatile approach to integrate heterogeneous datasets


% <img src="https://github.com/lauzingaretti/LinkHD/blob/master/linkhd.png" height="100" width="200">

## Laura Zingaretti & Yuliaxis Ramayo Caldas

#### m.lau.zingaretti@gmail.com

If you find our package useful, please cite:

Zingaretti,LM, Renand G, Morgavi DP, Ramayo-Caldas, Y (2019) LinkHD: a versatile framework to explore and integrate heterogeneous microbial comunnities. submitted. 

LinkHD is a general R software to integrate heterogeneous dataset focusing on micribial communities. LinkHD combines multivariate techniques to perform data integration with cluster and variable selection.
The method also allows us to study the relashionships between observations and features and to obtain **enrichment taxa analysis**. 

## Installation:

Clone the repository or 

devtools::install_github(repo="lauzingaretti/LinkHD")

## Usage

library(linkHD)

LinkHD capabilities were demostrated analizing public datasets from TARA Ocean expedition (https://oceans.taraexpeditions.org/en/m/about-tara/les-expeditions/tara-oceans/) and datasets from rumen metataxonomic communities (including bactera, archaea and protozoa data). Data can be loaded by the next command: 

data(Rumynotipes)

data(Tataoceans)

More examples and explanation of methods are available at:  https://lauzingaretti.github.io/LinkHD/
