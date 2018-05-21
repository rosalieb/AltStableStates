# Load datasets
# You can prepare your datasets in this script
# You will need a dataframe with your responses variable and your drivers
# If you work from multivariate data you need to transform them into univariate data using PCA or DCA for instance

# You can use our script with our Example data, but we also include below the code we used to clean it
# This code is in the 'prepare_datasets.R' script

## Example data ##
response <- read.delim(paste(getwd(),"/Input/Data/Bruel_etal_2018_response_example.txt",sep=""))
driver <- read.delim(paste(getwd(),"/Input/Data/Bruel_etal_2018_driver_example.txt",sep=""))
