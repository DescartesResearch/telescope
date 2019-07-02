# Telescope 
Telescope is a hybrid multi-step-ahead forecasting approach based on time series decomposition.

## Installation
This package can be installed in R by using the following commands:

`install.packages("devtools")` <br />
`devtools::install_github("DescartesResearch/telescope")` <br />

For unknown reasons, install_gitub does not work under all Windows versions. Therefore the package can alternatively be installed in R with the following commands:

`install.packages("remotes")` <br />
`remotes::install_url(url="https://github.com/DescartesResearch/telescope/archive/master.zip", INSTALL_opt= "--no-multiarch")`

## Example
`library(telescope)` <br />
`forecast <- telescope.forecast(taylor, horizon = 1000)`

For more information on this forecasting method, please visit our [homepage](http://descartes.tools/telescope).

![alt text](https://se.informatik.uni-wuerzburg.de/fileadmin/_processed_/7/3/csm_Telescope_982b20e78b.png "Telescope")
