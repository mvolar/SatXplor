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

SatXplor is configured to run on any Linux system. However there are a couple of dependencies that need to be installed before successful running of the pipeline.

0. Linux distribution
1. NCBI Blast
2. MAFFT multiple alignment 
3. R and some R package dependacies: 
   1. cmake
   2. libxml2-dev

## Installation

Since SatXplor is a compilation of over 10 individual scripts, and a Rust binary, the most elegant way of distributing and running EuSatXplor is by directly clonning the repository and installing the dependancies. 


1. Clone the repository:

    ```bash
    git clone https://github.com/mvolar/SatXplor.git
    ```

2. Navigate to the project directory:

    ```bash
    cd SatXplor
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

    Since R package manager in linux requires compilation of many packages, the installation time for all can take up to 20 minutes. Thus it is best to create a conda/mamba virtual environment and use the precompiled R packages for your Linux distirbution:
    
    ```
    mamba create -n myenv r-base=4.3.3 -c conda-forge -y 

    mamba activate myenv

    mamba install -c conda-forge r-biocmanager r-ggplot2 r-data.table r-dplyr r-umap r-stringr r-factominer r-ape r-optparse r-htmlwidgets r-igraph r-networkd3 r-circlize r-pheatmap r-scales bioconda::bioconductor-biostrings bioconda::bioconductor-complexheatmap  -y 
    ```


6. (Optional) Install both MAFFT and NCBI-BLAST through mamba:

    ```
    mamba activate myenv
    mamba install -c conda-forge -c bioconda mafft blast
    ```

## Docker

Alternative to the normal installation SatXplor also comes with a docker container:

1. Pull the containter
```
docker pull mvolaric/satxplor
```
2. Run the interactive shell and mount your data directory to `/mnt/data`
```
docker run -it -v path/to/your/data_folder:/mnt/data satxplor
```
3. In the interactive shell setup your desired `run_config.json` by using the provided helper script `setup_docker_run.py`:

```
python satxplor/setup_docker_run.py --input_genome_path genome.fasta \
--sat_raw satellites.fasta \
--final_results_dir output_folder
```

4. Run the `controller.py` from the interactive shell. After the run the results should be visible in `/path/to/your/data_folder/output_folder` directory.

```
python satxplor/controller.py
```

## Usage

Running SatXplor is simple, you just edit the `run_config.json` file:

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

And just run the main `controller.py` file which then runs and outputs the results in the `FINAL_RESULTS_DIR` and `SatXplor/results`.

```
python satxplor/controller.py
```

## Tests

To check if everything is installed correctly run the following scripts. The script searches for external dependencies (MAFFT, BLAST) as well as tests for all R import libraries.
``` bash
python tests/tests.py
```

SatXplor also ships whith a small testing sample to see if everything runs normally:

```
python satxplor/run_all_tests.py
```
which runs on the `testing_data/testing_data.tar.gz` files and produces the normal output of running EuSatXplor.

