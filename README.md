# Pathogen Island Detection Toolkit

## Description

This toolkit is designed for genomic sequence analysis and visualization, specifically for detecting pathogen islands in novel *Yersinia* strains. Developed as part of a university project, it aims to identify and analyze specific genes of interest in these bacterial strains.

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

These genes were referenced from a scientific article. [Include citation if available.]

## Features

This all-in-one solution operates from a single file, utilizing a collection of Python functions stored in a separate script. No coding experience is necessary – just follow the straightforward installation process and instructions to unlock a range of tasks related to gene analysis effortlessly. When executing the program, you'll be greeted by an intuitive menu presenting the following options: 

- **Gene Detection**: Determine if four specific pathogenic genes are present in the new *Yersinia* strains using BLASTP. For a detailed description of the output statistics, see [Output_Description.md](Output_Description.md).
- **Protein Search**: Analyze and determine the presence of any protein is in any strain. For a detailed description of the output statistics, see [Output_Description.md](Output_Description.md).
- **Sequence Retrieval**: Retrieve specific protein sequences by ID from a FASTA files.

## Installation

If using ubuntu, open ubunto files folder.

To install the required dependencies manually:
    ```bash
    pip install -r requirements.txt
    ```

Automated installation:
1. Clone the repository:
    ```bash
    git clone https://github.com/sapir-mardan/pathogen-genomic-analysis-toolkit.git
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



  
### Important Notes - MUST READ 

How to use FASTA files in the program: 
1. Always store protein sequence FASTA files in the designated "data" folder insidie "program_group6" folder. 
2. Ensure that your FASTA files have proper extensions such as ".fasta" or ".faa". 
3. When entering the name of the file to the program, include the extension.
4. Avoid using the same filename for multiple files to prevent conflicts.

Remove results folders after each use:  
  
When you use option 1 or option 2, your results would be created in a folder inside "code" folder.   
You must move the results folder with all its content from "code" folder after everytime you get results.  
If you dont, the next time you try and have results, they will be added to previous results and/or rewrtie the previous results.

# ACKNOWLEDGMENTS #  
  
This project was created by the amazing contributers:  
  
Sapir developed the python and bash script parts regarding the structure and running of the program. She also wrote the manual on how to properly install and run the program (README). Her expertise ensures a user-friendly experience, guiding users through a smooth program execution.  

Angela developed the code for running the scripts using programs such as blast. She also ..... (results for analysis)  

Both of them made the code compatible for both mac and windows users.  
    
  

