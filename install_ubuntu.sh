#!/bin/bash

# Verify current directory is pathogen-genomic-analysis-toolkit
if [ "$(basename "$PWD")" != "pathogen-genomic-analysis-toolkit" ]; then
    echo "Error: Please navigate to the 'pathogen-genomic-analysis-toolkit' directory before running this script."
    exit 1
fi

# Step 1: Install BLAST
echo "Installing BLAST..."
BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz"
curl -O $BLAST_URL
tar -xzf ncbi-blast-2.15.0+-x64-linux.tar.gz
mv ncbi-blast-2.15.0+ blast
rm ncbi-blast-2.15.0+-x64-linux.tar.gz

# Step 2: Install EDirect
echo "Installing EDirect..."
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
mv $HOME/edirect blast

# Step 3: Install required Python libraries
echo "Installing Python libraries..."
pip install biopython pandas

# Step 4: Verify installations
echo "Verifying installations..."
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

echo "Installation completed successfully."
