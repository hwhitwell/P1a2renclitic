## Synopsis

Biomarker networks based topological analysis of parenclitic networks. The networks are constructed by identifying pairs of variables
that deviate from a general population. This is achieved through the use of 2D kernal density estimation. The networks must be built using
continuous data, but can incorporate categorical data as well. Paper submitted (Whitwell 2017). A single network in build per individual
and predictive models are built based on topological features of the networks.

## Installation

The code can be run in R. "HW_Categorical_Parencltic.R" contains the code for the main program. In here, variables can be set and the file loaded. Detailed instructions are provided in the annotation. 
The folder "Functions" should be present in the same directry as HW_Categorical_Parenclitic. This contains the functions that "HW_Categorical_Parenclitic" calls. Data should be in a .csv format and
put into "Data" folder.

The script was built and tested using R3.4.0 in RStudio 1.0.143

Data: This should be in a .csv format and be layed out as follows:

first columns -> General infomation, case/control categories etc.
middle columns -> Continuous variables
end columns -> categorical variables if present.

There should be no missing data. If missing data is present, these rows/columns should be removed or the data imputed.

## Tests

An example data set (synthetic data) and results are provided. For infomation regarding the data, see (Whitwell 2017 - submitted).

## Contributors

Dr Harry J Whitwell
Dr Oleg Blyuss
Dr John F Timms
Prof Alexey Zaikin

## License

The code is open source, however we request that all use and reproduction should be referenced to Whitwell 2017 (submitted).