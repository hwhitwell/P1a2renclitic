## Synopsis

Biomarker modelling based on topological analysis of parenclitic networks. The networks are constructed by identifying pairs of variables
that deviate from a general population. This is achieved through the use of 2D kernal density estimation. The networks must be built using
continuous data, but can incorporate categorical data as well. Paper submitted (Whitwell 2017). A single network is built per individual
and predictive models are constructed based on topological features of the networks.

## Installation

The code is implemented in R. "HW_Categorical_Parencltic.R" contains the code for the main program. In here, variables can be set and the file loaded. Detailed instructions are provided in the annotation. 
The folder "Functions" should be present in the same directory as HW_Categorical_Parenclitic. This contains the functions that "HW_Categorical_Parenclitic" calls. Data should be in a .csv format and
put into "Data" folder.

The libraries that are required are listed in Functions/functions.R.

The script was built and tested using R3.4.0 in RStudio 1.0.143 on Windows 10.

Data: This should be in a .csv format and be layed out as follows:

first columns -> General infomation, case/control categories etc.
middle columns -> Continuous variables
end columns -> categorical variables if present.

There should be no missing data. If missing data is present, these rows/columns should be removed or the data imputed.

## Output files

A number of output files are generated.

Connections_XXX.csv = the number of connections each variable has for each individual.

Connections_XXX.png = a plot of the number of connections each variable has for cases and controls.

Contour_XXX_GridSize_YY.pdf = Networks for each individual.

IndexValues_XXX_GridSize_XX.csv = The topological index values for each index for each individual. These are used to build logistic regression models.

Results_XXX_GridSize_XX.pdf = XY plots of all combinations of pairs of topological indexes coded by case/control.

aucii.csv = Results of regression models. If multiple thresholds are run, the best result for each threshold will be appended to the file.

In aucii.csv, regression models are built based on topological indexes. The indexes in the models are indicated by a numerical code:

3 = Max degree

4 = Mean degree

5 = Network efficiency

6 = Max betweeness

7 = Mean betweeness

8 = Max closeness

9 = Mean closeness

10 = Max Google PageRank score

11 = Mean Google PageRank score

12 = Max Kleinberg's centrality score

13 = Mean Kleinberg's centrality score

14 = Total number of connections

15 = Max distance

16 = Mean distance

17 = Total number of connections to best marker

18 = Total number of connections to second best marker

19 = Total number of connections to third best marker

## Tests

An example data set (synthetic data) and results (20160615_Density) are provided. For infomation regarding the data, see (Whitwell 2017 - submitted).

## Contributors

Dr Harry J Whitwell

Dr Oleg Blyuss

Dr John F Timms

Prof Alexey Zaikin

## License

The code is open source, however we request that all use and reproduction should be referenced to Whitwell 2017 (submitted).