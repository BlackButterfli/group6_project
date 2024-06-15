# Pathogen Island Detection Toolkit

## Description

This toolkit is designed for genomic sequence analysis and visualization, specifically for detecting pathogen islands in novel *Yersinia* strains. Developed as part of a university project, it aims to identify and analyze specific genes of interest in these bacterial strains.

### Main Features:
The assigment asked to use BLAST spesifically in commend line. I have combined my experience with Python to create an easy-to-use interface that combines commend line and Python.  
  
_Bash Script Features:_  
**1. Retrieve Protein Sequences:** Automatically fetches protein sequences from the NCBI protein database using specified protein IDs.
**2. Merges multiple FASTA files.**
**3. Create BLAST Database from FASTA files.**:
**4. Perform BLAST Searches.**: 
**5. Clean Up:** Automatically removes intermediate files to keep the directory clean and organized.

Python Script :
The Python script integrates seamlessly with the bash script to run specific bioinformatics tasks.
Key functions that enhance the functionality of the pipeline:
**1.Use of Biopython Modules For Various Tasks (SeqIO, Entrez, NCBIXML):** 
    - Retrieve protein sequences from NCBI.
    - Parse and process BLAST results and save in a user-friendly format.
**2.Save and Adjust Detailed BLAST Results.**
**3.Sequence Retrieval from FASTA Files.** 


## Table of Contents

- [Project Background](#project-background)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

## Project Background

This toolkit was developed to support biologists in their research with data analyses and tools developments. They have recently isolated and sequenced a new strain of Yersinia bacteria and they suspect that 4 genes are involved in a new phenotype observed in their samples. The assignment involved preprocessing genomic data, performing sequence alignment, and visualizing the results to determine if pathogen islands are present in the given strains. 

The 4 genes were referenced from a scientific article and termed yRz1 yRz elyY holY:

Katharina Springer et al. "Activity of a Holin-Endolysin System in the Insecticidal Pathogenicity Island of Yersinia enterocolitica." *Journal of Bacteriology*, edited by Victor J. DiRita, Michigan State University. DOI: [10.1128/JB.00180-18](https://doi.org/10.1128/JB.00180-18)

## Features

This all-in-one solution operates from a single file, utilizing a collection of Python functions stored in a separate script. No coding experience is necessary – just follow the straightforward installation process and instructions to unlock a range of tasks related to gene analysis effortlessly. When executing the program, you'll be greeted by an intuitive menu presenting the following options: 

- **Gene Detection**: Determine if four specific pathogenic genes are present in the new *Yersinia* strains using BLASTP. For a detailed description of the output statistics, see [Output_Description.md]([Output_Description.md](https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit/blob/main/Output_Description.md)).
- **Protein Search**: Analyze and determine the presence of any protein is in any strain, and output the results in a format similar to 'Gene Detection' option.
- **Sequence Retrieval**: Retrieve specific protein sequences by ID from a FASTA files.

## Installation

**Automated Installation**
1. Clone the repository and change into the project directory:
    ```bash
    git clone https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit.git
    cd pathogen-genomic-analysis-toolkit
    ```
2. Change into the project directory:
    ```bash
    cd pathogen-genomic-analysis-toolkit
    ```
3. Run the installation script:
    ```bash
    ./install_mac.sh
    ./install_ubuntu.sh
    ```
4. If needed make script executable:
    ```bash
    chmod +x install_mac.sh
    chmod +x install_ubuntu.sh
    ```
    
**Manual Installation**
1. Clone the repository and change into the project directory:
    ```bash
    git clone https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit.git
    cd pathogen-genomic-analysis-toolkit
    ```

2. Install BLAST:
    ```bash
    BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-macosx.tar.gz"
    curl -O $BLAST_URL
    tar -xzf ncbi-blast-2.15.0+-x64-macosx.tar.gz
    # Move the extracted folder to the project directory and rename it to 'blast'
    mv ncbi-blast-2.15.0+ blast
    # Clean up by removing the tarball
    rm ncbi-blast-2.15.0+-x64-macosx.tar.gz
    ```

3. Install EDirect:
    ```bash
    sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
    echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile
    source $HOME/.bash_profile
    # Move the folder to blast folder directory
    mv $HOME/edirect blast
    ```

4. Install required Python libraries:
    ```bash
    pip install biopython pandas
    ```

5. Verify installations:
    ```bash
    cd blast
    if [ -d "bin" ] && [ -d "edirect" ]; then
        echo "BLAST and EDirect installed successfully."
    else
        echo "Error: BLAST or EDirect installation failed."
        exit 1
    fi

    cd ../code
    python3 -c "from Bio import Entrez; from Bio import SeqIO; from Bio.Blast import NCBIXML; import pandas as pd; import os; import sys"
    if [ $? -eq 0 ]; then
        echo "Python libraries installed successfully."
    else
        echo "Error: Python library installation failed."
        exit 1
    fi
    ```

If there are no errors, the installation is successful.

## Usage

1. Clone the repository (if not already done):
   ```bash
   git clone https://github.com/your-username/pathogen-island-detection-toolkit.git
   ```

2. Change into the project directory:
   ```bash
   cd pathogen-island-detection-toolkit
   ```

3. Run the program:
   ```bash
    ./code/run_program.sh
   ```

## Contribution 

Contributions are welcome!

**Acknowledgements:**    
**Sapir Mardan** developed the Python and Bash scripts for the program structure and execution. She also authored the README and the macOS manual, ensuring macOS compatibility.

**Angela Hsu** developed code for  Bash scripts insights on running BLAST and ensured Windows compatibility.

**Nicole Flores** and **Matthew Raj** provided valuable assistance throughout various sections of this project not presented here. Their help is greatly appreciated.

**University of Bristol** for their support and resources throughout this project.

## License

This project is licensed under the MIT License.

