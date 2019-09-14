# RT-qPCR_Plot

**To analyze qPCR results and generate bar plot accordingly**

## Raw data

To make it simple, <code>.csv</code> is the only valid format for data input. The sample data file is generated using **Applied Biosystems ViiA 7**. Should there be other systems being used, make sure the columns for ***experimental groups***, ***genes of interest***, and ***cycle threshold (ct) values*** are labeled as ***"Sample Name"***, ***"Target Name"***, and ***"CT"***, respectively. Name of ***control group*** as well as that for ***reference gene*** are expected to complete the following data analysis. Also, ***number of parallel ct values (experimental repeat)*** needs to be pointed out.

If the qPCR experiment is performed with an experimental repeat of 3, ct outliers should be identified and removed before analyzing the data for bar plot. Here, <font color='red'>**the most distant ct value apart from the middle one with a delta ct greater than 1 (2 folds)**</font> is set as the threshold for outlier identification.
