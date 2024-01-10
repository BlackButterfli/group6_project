# INTRODUCTION #

This program is runs from one file, utilizing a set of Python functions stored in a separate script to perform various tasks related to gene analysis. 
When executing the program, you have a main menu with 4 options:

Check for 4 Genes in New Yersinia (Option 1):
  * This option allows you to see how we extracted results (and extract yourself) about the four specific genes and if are present in new Yersinia.

Check Any Protein in Fasta File (Option 2):
  * You can use this option to check any protein of your own in any FASTA file within the "data" folder, and have similar results like in option 1.

Retrieve Sequence by ID (Option 3):
  * You can retrieve a sequence based on its ID from a specified FASTA file.

Exit (Option 'end'):
  * At any point you can type in "end" to go back to the main menu, or from the main menu to exit the program.


# INSTALLATION #

**For mac users:** (windows users skip)

### Step 1: Download Program Folder
  1. Download the folder "program_group6" and name it exactly "program_group6". Save where desired.
  2. Make sure you have 2 folders inside program_group6 folder: "code" and "data".
  3. Make sure you have 2 files inside code folder: "program_functions.py" and "run_program.sh".
  4. Make sure you have 2 files inside data folder: "NewYersinia1.faa" and "NewYersinia2.faa".

DO NOT CHANGE ANY NAMES OF FOLDER AND FILES, DO NOT MOVE ANY FOLDER OR FILE.

Personal note on this step: The reason the program is ordered specifically is to ensure both mac and windows can run it, and to make it a one click running program
(which means, everything you need is in the same code file, no need to maneuver between codes).
In fact, except from necessary installations, the user does not need to know or understand any coding. Only follow easy instructions,
Such as "what is your file name". The program will do the rest behind the scenes.


### Step 2: Install Python
  
  If you don't have Python installed, you can download it from the official Python website. 
  https://www.python.org/downloads/
  During installation, make sure to check the box that says "Add Python to PATH."

### Step 3: Install blast

=== explanation about blast

  1. Go to page: 
  https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
  2. Right click on "ncbi-blast-2.15.0+-x64-macosx.tar.gz". If needed choose "Save link as…" from the popup menu and save where desired.
  3. After finishing downloading, press on "ncbi-blast-2.15.0+-x64-macosx.tar" and a folder named "ncbi-blast-2.15.0+" will be extracted.
  4. Change the folder name to exactly "blast" and move it to the folder "program_group6".

### Step 4: Open Terminal
  In the Finder , open the /Applications/Utilities folder, then double-click Terminal.


### Step 5: Install EDirect

  To install the EDirect software, open a terminal window and execute one of the following two commands:

  sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

  sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
  
  Please allow 20 minutes for it to donwload. This will download a number of scripts and several precompiled programs into an "edirect" folder in the user's home directory. 
  It may then print an additional command for updating the PATH environment variable in the user's configuration file. The editing instructions will look something like:

  echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile

  You should get a message similar to: "Entrez Direct has been successfully downloaded and installed."

### Step 6: Move EDirect To Program File

  1. Open terminal.
  2. Type in:
     cd
  No messege should apear.
  3. Now you will copy to path to where EDirect is saved.
     Type in:
     find . -type d -name "edirect" 2>/dev/null
  4. Copy and save for later the path you got from after the ./
     Example: If you got: ./edirect -Copy: edirect
  5. Now you will copy to path to your blast folder in program folder:
     Type in:
     find . -type d -name "blast" 2>/dev/null

     Notice your path will end with /program_group6/blast
  6. Copy and save for later the path you got from after the ./
     Example: If you got: ./Desktop/bioinformatics/group_project/program_group6/blast -Copy: Desktop/bioinformatics/group_project/program_group6/blast
  7. Type in:
     mv (path from 4) (path from 6)
     For example, for me:
     mv edirect Desktop/bioinformatics/group_project/program_group6/blast

### Step 7: Install Required Libraries
  1. Open terminal and type:
     cd
  2. Type:
     pip install biopython pandas

### Step 8: Verify Installations
  1. Go to blast folder inside program_group6 folder. Make sure you have inside a folder named "bin", with many files inside, and a folder named "edirect" (among others).
  
  2. a. Open terminal and type:
        python3
     b. Type each line, one by one, press enter after each line:

        from Bio import Entrez
        from Bio import SeqIO
        from Bio.Blast import NCBIXML
        import sys
        import pandas as pd
        import os

     If there are no error messages, the libraries are successfully installed.

# RUN PROGRAM #

1. Open terminal and type
   cd
2. Type in:
   find . -type d -name "code" 2>/dev/null

   Again, Copy and save for later the path you got from after the ./
   Notice your path will end with /program_group6/code
   
3. Type in:
   bash run_program.sh

You will now be welcomed in the main menu.

# USING THE PROGRAM (PIPELINE) #

Using the program is streight forward. 

Option 1
Results saved in foler insise code.
After each use, move/delete the folder since another run will add to text.

Option 2
Results saved in folder inside code.
After each use, move/delete the folder since another run will add to text and overwrite csv.

Option 3


### Important Notes - must read

Whenever you are asked, ensure that the FASTA files for protein sequences are stored in the "data" folder.
A fasta file must end with fasta or faa

# ACKNOWLEDGMENTS #



פ

    
  

