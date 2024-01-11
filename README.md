# INTRODUCTION #

This program runs from one file, utilizing a set of Python functions stored in a separate script to perform various tasks related to gene analysis. 
When executing the program, you have a main menu with 4 options:

**Check for 4 Genes in New Yersinia (Option 1):**
  * This option allows you to see how we extracted results (and extract yourself) about the four specific genes and if are present in new Yersinia.

**Check Any Protein in Fasta File (Option 2):** 
  * You can use this option to check any protein of your own in any FASTA file and have similar results like in option 1.

**Retrieve Sequence by ID (Option 3):**
  * You can retrieve a sequence based on its ID from a specified FASTA file.

**Exit (Option 'end'):**
  * At any point you can type in "end" to go back to the main menu, or if you are in the main menu - enter ׳end׳ to exit the program.


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

### Step 3: Install Blast  

=== explanation about blast ===

  1. Go to page:   
  https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/  
  2. Right click on "ncbi-blast-2.15.0+-x64-macosx.tar.gz". If needed choose "Save link as…" from the popup menu and save where desired.
  3. After finishing downloading, press on "ncbi-blast-2.15.0+-x64-macosx.tar" and a folder named "ncbi-blast-2.15.0+" will be extracted.
  4. Change the folder name to exactly "blast" and move it to the folder "program_group6".

### Step 4: Open Terminal  
  In the Finder , open the /Applications/Utilities folder, then double-click Terminal.
  In a terminal, you simply type commands and press Enter to interact with a computer's operating system.
  We will use it to run the program.


### Step 5: Install EDirect  

  To install the EDirect software, open terminal (step 4) type one of the following and press Enter:  

  sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"  

  sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"  
    
  Please allow 20 minutes for it to donwload.   
  It may then print an additional command asking for PATH environment variable in the user's configuration file. The editing instructions will look something like:    
   
  echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile    
  
  If so, press Enter.  

  Now you should get a message similar to: "Entrez Direct has been successfully downloaded and installed."  

### Step 6: Move EDirect To Program File  

  1. Open terminal (step 4).  
  2. Type in and press Enter:  
     cd  
  No messege should apear.  
  3. Now you will copy to path to where EDirect is saved.    
     Type in and press Enter:   
     find . -type d -name "edirect" 2>/dev/null  
  4. Copy and save for later the path you got from after the ./    
     Example: If you got: ./edirect -Copy: edirect  
  5. Now you will copy to path to your blast folder in program folder:    
     Type in and press Enter:    
     find . -type d -name "blast" 2>/dev/null    
    
     Notice your path will end with /program_group6/blast  
  6. Copy and save for later the path you got from after the ./  
     Example:  
     If you got: ./Desktop/bioinformatics/group_project/program_group6/blast  
     -Copy:        Desktop/bioinformatics/group_project/program_group6/blast    
  8. Type in and press Enter:  
     mv (path from 4) (path from 6)  
     For example, for me:  
     mv edirect Desktop/bioinformatics/group_project/program_group6/blast  
 
### Step 7: Install Required Libraries  
  1. Open terminal (step 4), type and press Enter:  
     cd  
  2. Type and press Enter:  
     pip install biopython pandas  

### Step 8: Verify Installations  
  1. Go to blast folder inside program_group6 folder. Make sure you have inside a folder named "bin", with many files inside, and a folder named "edirect" (among others).  
    
  2. a. Open terminal (step 4), type and press Enter:  
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

1. Open terminal (step 4), type and press Enter:
   cd  
2. Type in and press Enter:  
   find . -type d -name "code" 2>/dev/null  
  
   Again, Copy and save for later the path you got from after the ./  
   Notice your path will end with /program_group6/code  
     
3. Type in and press Enter:   
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
Mention that it is NOT case sensitive.

  
### Important Notes - must read  

How to use FASTA files in the program: 
1. Always store protein sequence FASTA files in the designated "data" folder insidie "program_group6" folder. 
2. Ensure that your FASTA files have proper extensions such as ".fasta" or ".faa". 
3. When entering the name of the file to the program, include the extension.
4. Avoid using the same filename for multiple files to prevent conflicts.

# ACKNOWLEDGMENTS #  
  
This project was created by the amazing contributers:  
  
Sapir developed the python and bash script parts regarding the structure and running of the program. She also wrote the manual on how to properly install and run the program (README). Her expertise ensures a user-friendly experience, guiding users through a smooth program execution.  

Angela developed the code for running the scripts using programs such as blast. She also ..... (results for analysis)  

Both of them made the code compatible for both mac and windows users.  
    
  

