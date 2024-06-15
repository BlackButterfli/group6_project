#!/bin/bash

# Download URL for blast based on OS
if [ "$(uname -s)" == "Darwin" ]; then
    URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-macosx.tar.gz"
    curl -L $URL -o blast.tar.gz #Download 
elif  [ "$(Linux -s)" == "Darwin" ]; then
    URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz"
    wget $URL -O blast.tar.gz #Download 
else
    echo "Unsupported OS: $(uname -s)"
    echo "Blast was not installed"
    exit 1
fi

# Step 1: Install BLAST
echo "Installing BLAST..."
tar -xzvf blast.tar.gz    #unzip
sudo mv ncbi-blast-2.15.0+/* "\$HOME/blast"
export PATH=$PATH:"$HOME/blast/bin" # Add the BLAST bin directory to home directory

rm -rf ncbi-blast-2.15.0+ blast.tar.gz # Remove files
# Verify installations
if which blastn >/dev/null 2>&1; then
    echo "BLAST is installed."
else
    echo "BLAST is not installed."
fi

# Step 2: Install EDirect
echo "Installing EDirect..."
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=$PATH:"$HOME/edirect"


# Verify installations
if which efetch >/dev/null 2>&1; then
    echo "edirect is installed."
else
    echo "edirect is not installed."
fi
