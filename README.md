# SatXplor

A satDNA analysis pipeline.

## Table of Contents

- [EuSatXplor](#eusatxplor)
  - [Table of Contents](#table-of-contents)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Tests](#tests)

## Prerequisites

EuSatXplor is configured to run on any Linux system. However there are a couple of dependencies that need to be installed before successful running of the pipeline.

0. Linux distribution
1. NCBI Blast
2. MAFFT multiple alignment 
3. R and some R package dependacies: 
   1. cmake
   2. libxml2-dev

## Installation

Since EuSatXplor is a compilation of over 10 individual scripts, and a  Rust binary, the most elegant way of distributing and running EuSatXplor is by directly clonning the repository.

1. Clone the repository:

    ```bash
    git clone https://github.com/mvolar/EuSatXplor.git
    ```

2. Navigate to the project directory:

    ```bash
    cd EuSatXplor
    ```

3. Create and activate a virtual environment (optional but HIGHLY recommended):

    ```bash
    python -m venv venv
    source venv/bin/activate  
    ```

4. Install python dependencies:

    ```bash
    pip install -r requirements.txt
    ```

5. Install R dependencies:

    Since R package manager is a bit tricky on linux go the R source files folder and open up your R session
    ```
    cd /eusatxplor/r/
    R
    ```
    Then run:

    ```{r}
    source("install_packages.R")
    ```
    An R session installer will guide you. If there are errors during installation, most likely explanation is you are missing some developmental libraries of C/C++ (like cmake, libxml2), all of which R will let you know that you are missing (if you install R )will then have to be installed by:

    ```
    sudo apt install "missing_lib"
    ```
    After installing simply rerun the `R/source()` command and the installation of R packages will continue

## Usage

Running EuSatXplor is simple, you just edit the `run_config.json` file:

```json
{
    "INPUT_GENOME_PATH": "path/to/your/genome",
    "GENOME_PATH": "./input.fasta", #this is the location of the temporary copy
	"SAT_RAW": "path/to/your/sats",
    "SAT_FASTA_PATH": "./sats.fasta", #this is the location of the temporary copy
    "FINAL_RESULTS_DIR": "/path/to/final/results/dir",
    "OVERWRITE": true
}
```

And just run the main `controller.py` file which then runs and outputs the results in the `FINAL_RESULTS_DIR` and `EuSatXplor/results`.


## Tests

To check if everything is installed correctly run the following scripts. The script searches for external dependencies (MAFFT, BLAST) as well as tests for all R import libraries.
``` bash
python tests/tests.py
```

EuSatXplor also ships whith a small testing sample to see if everything runs normally:

```
python eusatxplor/run_all_tests.py
```
which runs on the `testing_data/testing_data.tar.gz` files and produces the normal output of running EuSatXplor.

