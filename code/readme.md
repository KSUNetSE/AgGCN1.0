# Gene Co-expression Network Generation
This is a Python package for the generation of Gene Co-expression Networks.

# Requirements
Python environment 3.6.0+.

Install [Anaconda](https://www.anaconda.com/) (preferred) or [Python](https://www.python.org/downloads/). Check the box "add python command to path" when installing the package.

When the Python environment is installed, the following required packages can be installed using the `pip` command from the terminal:

`pip install numpy csv scipy matplotlib`

activate the base environment in the terminal:

`activate base`

# Inputs
The input co-expression dataset must be saved in .csv format, rows represent nodes and columns represent experimental conditions (1st column: gene names, 1st row: conditions). The dataset may contain missing values.

Run the following command to show the input parameters:

`python GCN.py --help`

Input parameters: 

`--data` 

path to the co-expression data matrix, e.g., `--data anopheles_expression_data.csv`

default: anopheles.csv

`--zero_removed` 

True: 0s have been removed from the input dataset or do not want to remove 0s from the input dataset. 

False: remove all 0s from the input dataset. 

default: True

`--zscored`

True: Each column has been z-score normalized or the user do not want to normalize the data.

False: z-score normalize each column.

default: True

`--rescaled`

True: each value has been log2 rescaled. This is used for RNA-seq data which are distributed in wide ranges.

False: for any value `x` in the input dataset, it will be rescaled with `log2(x)`.

default: True

`--dropna`

type: interger, keep only the conditions with at least `x` genes.

default: 200

`--save_data`

True: save the processed data as: 'z-scored xxx.csv', where 'xxx' represents the name of input dataset.

False: the processed data will not be saved.

default: True

`--curve_params`

True: use pre-determined parameters for the sliding threshold curve.

False: the parameters of the sliding threhold curve will be fitted from the PCCs.

default: False

`--alpha`

input pre-determined `alpha` if `curve_params` is True.

`--eta`

input pre-determined `eta` if `curve_params` is True.

`--lam`

input pre-determined `lam` if `curve_params` is True.

`--beta`

input pre-determined `beta` if `curve_params` is True.

`--pcc`

save the PCC matrix, e.g., `--pcc PCC_matrix.csv`.

default: void

`--paired_elements`

save a matrix which contains the number of paired elements between every gene pair. e.g., `--paired_elements paired_elements_matrix.csv`.

default: void

`--bin_size`

 binning the PCCs into different intervals per paired element, for example:
 
|4, 5, ... , 10| 11, 12, ... , 20 | 21, 22, ... , 30 | 
|----------------|-------------------------------|-----------------------------|

 default: 10
 
 `--cutoff`
 
 choose the top fraction of gene pairs with the highest PCCs as the edges of the network, e.g., `--cutoff 0.005`
 
 default: 0.005
 
 `--thres_curve`
 
 save the fitted sliding threshold curve, e.g., `--threshold_curve fitted_threhold_curve.png`.
 
 default: threshold_curve.png
 

# Default outputs:

`z-scored anopheles.csv`

`PCC_matrix.csv`

`paired_elements.csv`

`edgelist.csv`

`threshold_curve.png`

# Examples

## Example 1

We generate a gene co-expression network for the `anopheles.csv` dataset through the following command:

`python GCN.py --data anopheles.csv --zero_removed True --zscored True --rescaled True --save_data True --dropna 200 --pcc anopheles_PCC_matrix.csv --paired_elements anopheles_paired_elements.csv --edgelist anopheles_edgelist.csv --thres_curve anopheles_threshold_curve.png`

outputs:

`z-scored anopheles.csv`

`anopheles_PCC_matrix.csv`

`anopheles_paired_elements.csv`

`anopheles_edgelist.csv`

`anopheles_threshold_curve.png`

## Example 2

We generate a gene co-expression network for the `aedes.csv` dataset through the following command:

`python GCN.py --data aedes.csv --zero_removed False --zscored False --rescaled False --save_data True --dropna 200 --pcc aedes_PCC_matrix.csv --paired_elements aedes_paired_elements.csv --edgelist aedes_edgelist.csv --thres_curve aedes_threshold_curve.png`

outputs:

`z-scored aedes.csv`

`aedes_PCC_matrix.csv`

`aedes_paired_elements.csv`

`anopheles_edgelist.csv`

`aedes_threshold_curve.png`


# Citation

A previous version of the code can be found in the [paper](https://doi.org/10.1186/s12859-022-04697-9):

Kuang, J., Buchon, N., Michel, K. _et al._ A global Anopheles  gambiaeAnopheles gambiae gene co-expression network constructed from hundreds of experimental conditions with missing values. _BMC Bioinformatics_  23, 170 (2022). https://doi.org/10.1186/s12859-022-04697-9

