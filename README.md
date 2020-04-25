# Telescope 
Telescope is a hybrid multi-step-ahead forecasting approach based on time series decomposition. 

Details of the Telescope approach can be found at [1,2]. The details of the recommendation approach can be found at [3].

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


[1] Bauer, A., Züfle, M., Herbst, N., Kounev, S. & Curtef, V. (2020). [Telescope: An Automatic Feature Extraction and Transformation Approach for Time Series Forecasting on a Level-Playing Field](https://www.bibsonomy.org/documents/8efb3f8c174e0904cce5bdaadb3e6160/andre.bauer/BaZuHeKoCu-ICDE-Telescope.pdf). Proceedings of the 36th International Conference on Data Engineering (ICDE) (p./pp. 1902-1905).

[2] Bauer, A., Züfle, M., Herbst, N., Zehe, A., Hotho, A. & Kounev, S. (2020). [Time Series Forecasting for Self-Aware Systems](https://doi.org/10.1109/JPROC.2020.2983857). Proceedings of the IEEE. 

[3] Bauer, A., Züfle, M., Grohmann, J., Schmitt, N., Herbst, N. & Kounev, S. (2020). [An Automated Forecasting Framework based on Method Recommendation for Seasonal Time Series](https://www.bibsonomy.org/documents/27a54ab734579fafad2f3bc44eea8daf/andre.bauer/BaZuGrScHeKo-ICPE20-Seasonal-Forecast.pdf). Proceedings of the 2020 ACM/SPEC International Conference on Performance Engineering (p./pp. 48-55), April, New York, NY, USA: ACM.
