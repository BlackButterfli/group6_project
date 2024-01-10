#!/bin/bash

#For mac users we provide the script with alternative paths where the blast end edirect programs are saved. The script will automatically look for them, and that way the code is written the same for both mac and obuntu users.
export PATH=$PATH:../blast/bin
export PATH=$PATH:../blast/edirect

#############   OPTION 1 START  #######################

option_1()
{

##Checks NewYersinia Fasta Files against 4 Proteins in question

#Create a path and folder for results:
mkdir option1_results
pathway="option1_results/results.csv"

# Retrieves protein seuquences using protein IDS to create database
esearch -db protein -query "CAI77377.1 OR CAI77378.1 OR CAI77379.1 OR SOH98281.1" | efetch -format fasta > ../data/proteins.fasta

# Combines the two NewYersinia files into 1 fasta file
cat ../data/NewYersinia1.faa ../data/NewYersinia2.faa > ../data/CombinedYersinia.faa

# Makes database to perform blast
makeblastdb -in ../data/CombinedYersinia.faa -dbtype prot -out NewYersiniadb

#makes CSV file with output file destination
blastp -query ../data/proteins.fasta -db NewYersiniadb -out "$pathway" -evalue 0.001 -outfmt "10 qseqid sseqid pident length mismatch qstart qend evalue"

# Call the adjust csv function from python script
python3 program_functions.py adjust_csv "$pathway"

##############################
# Create XML and parse to txt
###############################
#Create a path for results:
name_txt="option1_results/results.txt"
name_xml=xml1

#makes XML file with output file destination
blastp -query ../data/proteins.fasta -db NewYersiniadb -out "$name_xml" -evalue 0.001 -outfmt 5

#Parse the XML and save as txt with functions from python script:
python3 program_functions.py print_results_blastp "$name_xml" "$name_txt"
 
# Deletes the extra files created to make database to keep directory clear
find . -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +

#Print the user the results:
echo
echo Your results are saved in folder "option1_results".
echo You will have 2 types of results: 
echo 1. A csv file with summary data from blastp
echo 2. A txt file with aligment information from blastp
echo
echo ---------------------------------------------------------------------------------------------------------
}
#End option 1


#############   OPTION 2 START  #######################

option_2() {

#########################
#check fasta file exists
##########################

#A loop that asks fasta file name from user and checks if exists:

while true; do
python3 program_functions.py option_2
read -p "Enter the name of your FASTA file, do not use quotations: " file_path
in_pathway="../data/$file_path"

#Option to exit to main menu
if [ "$file_path" == "end" ]; then
echo You chose to return to the main menu
return
fi

#Checks if fasta name exists (e.g., -e)
if [ -e "$in_pathway" ]; then
echo Your fasta file is: $file_path
break
else
echo "The file '$file_path' does not exist. Please enter the name of an existing file."
continue
fi
done

# Makes database to perform blast from the fasta
makeblastdb -in "$in_pathway" -dbtype prot -out NewTestdb


##################################################################
# retrive protein ID and save it as fasta file "proteins.fasta"
##################################################################
#A loop that checks if a protein id exists in the fasta provided:

while true; do
#Will print from py: "Enter the Protein IDs you want to search for. If multiple protein IDs separate with OR (i.e. SOH98281.1 OR CAI77379.1)":
python3 program_functions.py "Print_User"
read protein_ids

#an option for user to go back to the main menu
if [ "$protein_ids" == "end" ]; then
echo You chose to return to the main menu
return
fi

# Retrieves protein seuquences using protein IDS to create database
esearch -db protein -query "$protein_ids" | efetch -format fasta > proteins.fasta

# checks whether the file temp.fasta exists and has a size greater than zero:
if [ -s proteins.fasta ]; then
break

else
echo "Invalid protein IDs. Please try again:"
fi
done
#End of the loop #

###########################
# CSV
##########################
#Create folder and path for results:
mkdir option2_results
out_pathway="option2_results/results.csv"

#makes CSV file with output file destination and adjust it with python function:
blastp -query proteins.fasta -db NewTestdb -out "$out_pathway" -evalue 0.001 -outfmt "10 qseqid sseqid pident length mismatch qstart qend evalue"
python3 program_functions.py adjust_csv "$out_pathway"

#############################
# THE XML
############################
#Create path for results:
path_xml="thefile.xml"
path_txt="option2_results/results2.txt"

#makes XML file with output file destination
blastp -query proteins.fasta -db NewTestdb -out $path_xml -evalue 0.001 -outfmt 5

#Parse the XML file and save as txt file
python3 program_functions.py print_results_blastp "$path_xml" "$path_txt"


# Deletes the extra files created to make database to keep directory clear
find . -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +

#Print user the results:
echo
echo Your results are saved in folder "option2_results".
echo You will have 2 types of results: 
echo 1. A csv file with summary data from blastp  
echo 2. A txt file with aligment information from blastp       
echo

}

#End option 2

#################################
#RUNNING THE BASH
#################################


while true; do
	python3 program_functions.py

	choise1=$(echo $?)
	
	if [ $choise1  == "1" ]; then
    		option_1
	fi

	if [ $choise1  == "2" ]; then
		option_2	
	fi

	if [ $choise1  == "3" ]; then
    		python3 program_functions.py option_3
	fi

	if [ $choise1  == "9" ]; then
    		break
	fi
done
find . -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +
echo
echo Goodbye!
exit 0
