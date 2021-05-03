# Robust-estimation-of-SARS-CoV-2-epidemic-in-US-counties

This repository contains codes and data for the paper [Robust estimation of SARS-CoV-2 epidemic in US counties](https://arxiv.org/pdf/2010.11514.pdf). It's separated into two parts:

1. Folder *Reproducing results in the paper* provides codes and data that could reproduce all the figures and tables in the paper.
2. Folder *Forecast in US counties* provides codes and data that can generate forcasts in US counties with more than 2 death cases as of September 20, 2020. The forecast start date and length can be easily adjusted by changing the variables "end_date" and "prediction_length" in the codes.

To run codes in the above two folders, please follow the steps below:

Step 1. Install [Rstudio](https://www.rstudio.com/) and update the [R software](https://www.r-project.org/) to the latest version (4.1.0 currently).

Step 2. Install all the R packages needed for the codes. Especially for the package "RobustGaSP", please make sure it's updated to the version 0.6.1 or above.

Step 3. Make sure the working directory of your Rstudio is under the master branch.

Step 4. Now you can run the R codes in these folders and reproduce the results in our paper!

As of the output, 

1. for codes in the folder *Reproducing results in the paper*, it will automatically produce EPS figures and save them to the folder *Results* under *Reproducing results in the paper*.
2. for codes in the folder *Forecast in US counties*, it will generate 51 RDatas containing estimations, forecasts and uncertainty quantifications of counties in 51 US states. All outputs will be saved in the folder *Results* under *Forecast in US counties*.

You are more than welcome to post issues or contact us if you have any questions!

Email: mengyang@pstat.ucsb.edu & hanmo@pstat.ucsb.edu
