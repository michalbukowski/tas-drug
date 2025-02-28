## A Jupyter notebook for co-occurrence analysis of toxin-antitoxin (TA) systems and drug resistance determinants (tas-drug)

Jupyter Python notebook for the [tblast-genomes](https://github.com/michalbukowski/tblastn-genomes) pipeline output analysis and visualisation. For a broader context see the following publication:

Bukowski M, Banasik M, Chlebicka K, Bednarczyk K, Bonar E, Sokołowska D, Żądło T, Dubin G, Władyka B. _Analysis of co-occurrence of type II toxin–antitoxin systems and antibiotic resistance determinants in Staphylococcus aureus_. mSystems 0:e00957-24.
https://doi.org/10.1128/msystems.00957-24

The analysis utilises CARD data (located in the input directory):

Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, Edalatmand A, Petkau A, Syed SA, Tsang KK, Baker SJ. _CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database_. Nucleic acids research. 2023 Jan 6;51(D1):D690-9.

The most current data set may be obtained from [card.mcmaster.ca](https://card.mcmaster.ca) (Download/Download CARD Data). The pipeline utilises the protein homology model sequences.

Subsequent steps of the analysis are described in detail in the `TAS-Drug.ipynb` notebook.

### 1. Environment
`Miniconda`/`Anaconda` installation is a prerequisite. The `tas-drug` environment must be created prior running Jupyter Lab and the `TAS-Drug.ipynb` notebook code. The repository has been tested on `Ubuntu 22.04` using `conda 24.11`.

In `envs` directory the `tas-drug.yml` file describes the conda environment suitable for the analysis. If you want to install the latest versions of the packages, remove all version designations from that file (`=X.X.*`). Run the following command from the repository directory to create the `tas-drug` environment:
```
conda env create -f envs/tas-drug.yml
```

### 2. The input data
In order to create the `input` directory and populate it with input files (deposited elsewhere) run the `get_input.sh` script:
```
bash get_input.sh
```
The `input` directory should be populated with files listed in the next section. The input files will contain the original input data that were used in the study by [Bukowski et al. (2025)](https://doi.org/10.1128/msystems.00957-24), including the `tblastn.tsv` file, which contains the output of the [tblast-genomes](https://github.com/michalbukowski/tblastn-genomes) pipeline.

### 3. Directory structure
```
tas-drug/
├── envs/
│   └── tas-drug.yml
├── lib/
│   ├── assembly.py
│   └── functions.py
├── input/
│   ├── genomes_list.txt
│   ├── card.faa
│   ├── antitoxins.faa
│   ├── toxins.faa
│   └── tlastn.tsv
├── output/
├── get_input.sh
└── TAS-Drug.ipynb
```

In the working directory, you will find the `envs` and `input` directories as well as the `get_input.sh` script, all of which have been described in the previous section. The `TAS-Drug.ipynb` notebook is the main engine for performing the analysis. The `lib` directory contains two files with the `Assembly`, `Locus` and `Gene` class definitions as well as definitions of 8 functions, all utilised by the `TAS-Drug.ipynb` notebook. The `output` directory will be automatically created when the analysis is run. The details of the analysis is described step by step in the `TAS-Drug.ipynb` notebook and more detailed information can be found by inspecting the functions docstrings.

For the reference, the `TAS-Drug.ipynb` contains complete original output of the analysis obtained for the original input data.

### 4. Running the analysis
To run the analysis activate the `tas-drug` environment and run Jupyter Lab in the repository directory:
```
conda activate tas-drug
jupyter lab
```
In the Jupyter Lab app open open the `TAS-Drug.ipynb` notebook and run the code. Make sure you are using a Python kernel.
