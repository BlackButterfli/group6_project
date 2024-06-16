# Pathogen Island Detection Toolkit

## Description

A toolkit designed for genomic sequence analysis, specifically for detecting pathogen islands in novel *Yersinia* strains. Developed as part of a university project, it aims to identify and analyze specific genes of interest in these bacterial strains.

### Main Features:
The assigment asked to use BLAST spesifically in commend line. I have combined my experience with Python to create an easy-to-use interface, with Python script that integrates seamlessly with the bash script to run specific bioinformatics tasks and enhance the functionality of the pipeline.
  
**Key Bash Script Features:**
  - Create BLAST database from FASTA files and perform BLAST searches.
  - Retrieve protein sequences from NCBI protein database.
  - Automatically removes intermediate files to keep the directory clean and organized.
  
**Key Python Script Features:**  
  - Use of Biopython Modules For Various Tasks (SeqIO, Entrez, NCBIXML).
  - Retrieve protein sequences from NCBI.
  - Sequence Retrieval from FASTA Files. 
  - Parse and process BLAST results and save output in a user-friendly format.

## Table of Contents

- [Project Background](#project-background)
- [Features](#features)
- [Installation & Requirements](#installation-and-requirements)
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

<img width="738" alt="image" src="https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit/assets/92859243/d529deb2-385e-4ddc-86be-205d9a642d78">


- **Gene Detection**: Determine if four specific pathogenic genes are present in the new *Yersinia* strains using BLASTP. For a detailed description of the output statistics, see [Output_Description.md]([Output_Description.md](https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit/blob/main/Output_Description.md)).
- **Protein Search**: Analyze and determine the presence of any protein is in any strain, and output the results in a format similar to 'Gene Detection' option.
- **Sequence Retrieval**: Retrieve specific protein sequences by ID from a FASTA files.

## Installation and Requirements

**Requirements**
BLAST, edirect, python3 and the following libraries:
Biopython (Bio)
sys
pandas
os
datetime

**Installation**
1. Clone the repository and change into the project directory:
    ```bash
    git clone https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit.git
    cd pathogen-genomic-analysis-toolkit
    ```
2. Change into the project directory:
    ```bash
    cd pathogen-genomic-analysis-toolkit
    ```
3. Run the installation script to download blast and edirect to your computer (linux/mac):
   ```bash
    ./install_blast_edirect.sh
    ```
   If needed make script executable:
    ```bash
    chmod +x install_blast_edirect
    ```
4. Download required libraries:
   ```bash
   pip install -r requirements.txt
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
**Sapir Mardan** developed the Python and Bash script for the program structure and execution.

**Angela Hsu** developed code for Bash script, specifically for running BLAST and ensured Windows compatibility.

**Nicole Flores** and **Matthew Raj** provided valuable assistance throughout various sections of this project not presented here. Their help is greatly appreciated.

**University of Bristol** for their support and resources throughout this project.

## License

This project is licensed under the MIT License.

