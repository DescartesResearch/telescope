# Telescope 
Telescope is a machine learningbased forecasting approach that automatically retrieves relevant information from a given time series.. 

![alt text](https://se.informatik.uni-wuerzburg.de/fileadmin/_processed_/7/3/csm_Telescope_982b20e78b.png "Telescope")

Details of the Telescope approach can be found at [[1](https://www.bibsonomy.org/documents/8efb3f8c174e0904cce5bdaadb3e6160/andre.bauer/BaZuHeKoCu-ICDE-Telescope.pdf),[2](https://doi.org/10.1109/JPROC.2020.2983857)]. The details of the recommendation approach can be found at [[3](https://www.bibsonomy.org/documents/27a54ab734579fafad2f3bc44eea8daf/andre.bauer/BaZuGrScHeKo-ICPE20-Seasonal-Forecast.pdf)].

## Installation
This package can be installed in R by using the following commands:

```
install.packages("devtools")
devtools::install_github("DescartesResearch/telescope")
```

For unknown reasons, install_github does not work under all Windows versions. Therefore the package can alternatively be installed in R with the following commands:

```
install.packages("remotes")
remotes::install_url(url="https://github.com/DescartesResearch/telescope/archive/master.zip", INSTALL_opt= "--no-multiarch")
```

## Example without Recommendition system
```
library(telescope)
forecast <- telescope.forecast(forecast::taylor, horizon = 1000)
```

## Example with Recommendition system
```
install.packages('Mcomp')
library(Mcomp)
library(telescope)

ts.list <- list()
for(i in 1:length(M3)){
  ts.list[i] <- list(ts(c(M3[[i]]$x,M3[[i]]$xx),frequency = frequency(M3[[i]]$x)))` <br />
}

model <- telescope.trainrecommender(ts.list)
telescope.forecast(forecast::taylor, horizon = 1000, rec_model = model)
```

## Further Information and References
For more information on this forecasting method, please visit our [homepage](http://descartes.tools/telescope).

[1] Bauer, A., Züfle, M., Herbst, N., Kounev, S. & Curtef, V. (2020). [Telescope: An Automatic Feature Extraction and Transformation Approach for Time Series Forecasting on a Level-Playing Field](https://www.bibsonomy.org/documents/8efb3f8c174e0904cce5bdaadb3e6160/andre.bauer/BaZuHeKoCu-ICDE-Telescope.pdf). Proceedings of the 36th International Conference on Data Engineering (ICDE) (p./pp. 1902-1905).

[2] Bauer, A., Züfle, M., Herbst, N., Zehe, A., Hotho, A. & Kounev, S. (2020). [Time Series Forecasting for Self-Aware Systems](https://doi.org/10.1109/JPROC.2020.2983857). Proceedings of the IEEE. 

[3] Bauer, A., Züfle, M., Grohmann, J., Schmitt, N., Herbst, N. & Kounev, S. (2020). [An Automated Forecasting Framework based on Method Recommendation for Seasonal Time Series](https://www.bibsonomy.org/documents/27a54ab734579fafad2f3bc44eea8daf/andre.bauer/BaZuGrScHeKo-ICPE20-Seasonal-Forecast.pdf). Proceedings of the 2020 ACM/SPEC International Conference on Performance Engineering (p./pp. 48-55), April, New York, NY, USA: ACM.
