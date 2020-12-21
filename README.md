This repository includes the codes (R and Stan) to apply multilevel regression and poststratification (MRP) to COVID test results, as a supplement to the paper "A Proxy Metric for Community Coronavirus Tracking", by Len Covello , Andrew Gelman , Yajuan Si , and Siquan Wang.


The Rmd file "hospital monitoring.Rmd" includes the R codes for data importing, cleaning, modeling, and generating graphics.


The Stan file "MRP.stan" has the Stan codes to apply MRP in the adjustment of selection bias.


The "testdata.csv" provides a format template of test results. The original data are confidential and not for public release. We need the individual-level test results and dates, and patients' characteristics (age, gender, race/ethinicty and county information).

The remaining datasets from the public census and state websites. We cannot release the data of the hospital population used for poststratification, which provide the number of patients by the crosstabulation of age, gender, race/ethinicty and couny. 
