#!/bin/bash


#############   OPTION 1  #######################

option_1()
{
	
##Checks NewYersinia Fasta Files against 4 Proteins in question

#Create a path and folder for results:
mkdir option1_results
pathway="option1_results/results.csv"

# Retrieves protein seuquences using protein IDS to create database
esearch -db protein -query "CAI77377.1 OR CAI77378.1 OR CAI77379.1 OR SOH98281.1" | efetch -format fasta > data/proteins.fasta

# Combines the two NewYersinia files into 1 fasta file
cat data/NewYersinia1.faa data/NewYersinia2.faa > data/CombinedYersinia.faa

# Makes database to perform blast
makeblastdb -in data/CombinedYersinia.faa -dbtype prot -out NewYersiniadb

#makes CSV file with output file destination
blastp -query data/proteins.fasta -db NewYersiniadb -out "$pathway" -evalue 0.001 -outfmt "10 qseqid sseqid pident length mismatch qstart qend evalue"

# Call the adjust csv function from python script
python3 program_functions.py adjust_csv "$pathway"
rm option1_results/results.csv

#makes XML file with output file destination
blastp -query data/proteins.fasta -db NewYersiniadb -out "xml1" -evalue 0.001 -outfmt 5

#Parse the XML and save as txt with functions from python script:
python3 program_functions.py  save_results_blastp "xml1" "option1_results"

# Deletes the extra files created to make database to keep directory clear
find . -maxdepth 1 -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +

}
#End option 1


#############   OPTION 2 START  #######################

option_2() {

#A loop that asks fasta file name from user and checks if exists:
while true; do
python3 program_functions.py option_2
read -p "Enter the name of your FASTA file, do not use quotations: " file_path
in_pathway="data/$file_path"

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

#Create folder and path for results:
mkdir option2_results
out_pathway="option2_results/results.csv"

#makes CSV file with output file destination and adjust it with python function:
blastp -query proteins.fasta -db NewTestdb -out "$out_pathway" -evalue 0.001 -outfmt "10 qseqid sseqid pident length mismatch qstart qend evalue"
python3 program_functions.py adjust_csv "$out_pathway"
rm option2_results/results.csv

#makes XML file with output file destination
blastp -query proteins.fasta -db NewTestdb -out xml2 -evalue 0.001 -outfmt 5

#Parse the XML file and save as txt file
python3 program_functions.py save_results_blastp "xml2" "option2_results"

# Deletes the extra files created to make database to keep directory clear
find . -maxdepth 1 -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +
rm data/proteins.fasta

}

#RUNNING THE PROGRAM
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
find . -maxdepth 1 -type f ! \( -name '*.py' -o -name '*.sh' -o -name '*.csv' -o -name '*.txt' \) -exec rm {} +
echo 
echo Goodbye!
exit 0


