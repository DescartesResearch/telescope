# Telescope 
Telescope is a hybrid multi-step-ahead forecasting approach based on time series decomposition.

## Installation
This package can be installed by using the following commands:

`install.packages("devtools")` <br />
`devtools::install_github("marwinzuefle/telescope")` <br />
`library(telescope)`

## Example
`forecast <- telescope.forecast(AirPassengers, horizon = 10)`

For more information on this forecasting method, please visit our [homepage](http://descartes.tools/telescope).

![alt text](https://se.informatik.uni-wuerzburg.de/fileadmin/_processed_/7/3/csm_Telescope_982b20e78b.png "Telescope")
