from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIXML
import sys
import pandas as pd
import os

########################################################################
#Code for user menu:
########################################################################

def user_main_menu():  
    """
    This function gives the user the option to choose between different options for the program.
    If the user didnt type one of the option, it will ask it for another try.
    If the user choose "end", the bash script will terminate.
    The bash script gets an exit value (sys.exit()) which it interpret as an argument, and from which knows how to continue.
    """  
    while True: 
        print('\033[4m' + "\nWelcome to the main menu.\n" + '\033[0m')
        print("How would you like to continue?\n")
        print("To check if 4 genes are in new Yersinia enter 1")
        print("To check any protein of your own in any fasta file data enter 2")
        print("To retrive a sequence based on its ID from a FASTA file enter 3")
        print( "\033[95m" + "To exit enter 'end'\n" + "\033[0m")
        user_choose = input("Enter your choice:")
        print(f"You chose option: {user_choose}")
        
        if user_choose == "1":
            #option_1()
            sys.exit(1)

        elif user_choose == "2":
            sys.exit(2)

        elif user_choose == "3":
            sys.exit(3)
            
        elif user_choose.lower() == "end":
            #print("\nGoodbye!")
            sys.exit(9)
            
        else:
            print("\nTry again!\n")        
        

        

##########################################
# Option 1 + 2
##########################################

def adjust_csv(pathway):
    """
    This function takes a csv output for spcific genes, and make it better in terms of visuallity.
    """
    # Load the CSV file using the provided file_pathway
    df = pd.read_csv(pathway, header=None)
    # Add column headers
    column_headers = ["Protein ID", "Sequence ID", "% identity", "Alignment Length", "Mismatches", "Query Start", "Query End", "E-Value"]
    df.columns = column_headers
    # Changes protein ID's to Protein names
    protein_names = {
        "CAI77377.1": "holY",
        "CAI77378.1": "elyY",
        "CAI77379.1": "yRz",
        "SOH98281.1": "yRz1"
    }
    df.replace(protein_names, inplace=True)
    # Save the modified DataFrame back to CSV using the same file path
    df.to_csv(pathway, index=False)

protein_names = {
    "CAI77377.1": "holY",
    "CAI77378.1": "elyY",
    "CAI77379.1": "yRz",
    "SOH98281.1": "yRz1"
}    

def print_results_blastp(output_file, txt_name):
    """
    This function parses BLASTP XML output file, prints alignment details, and saves results to a text file.

    Parameters:
    - output_file (str): Path to the BLASTP XML output file.
    - txt_name (str): Name of the text file to save the results.
    """
    
    # parse the XML file using NCBIXML module
    with open(output_file, "r") as results:
        blast_records = NCBIXML.parse(results)
        
        # print for the user about the results being printed and saved
        print("Results printed here are also saved as a txt file:")
        
        # iterates over records (we only have one so we can use read - I need to ask you what you prefer)
        for record in blast_records:
            # iterates over alignments in the each record
            for alignment in record.alignments:
                print(f"\nAlignment: {alignment.title}")
                if record.query.split(" ")[0] in protein_names:
                    print(f"ID {record.query.split(' ')[0]} is: '{protein_names[record.query.split(' ')[0]]}'")
                print(f"Record ID: {record.query}")
                #iterates over high-scoring pairs (HSPs) in the alignment
                for hsp in alignment.hsps:
                    #print detailed information about each HSP
                    print(f"  E-value:   {hsp.expect}")
                    print(f"  Bit Score: {hsp.bits}")
                    print(f"  Length: {alignment.length}")
                    print(f"  Number of gaps: {hsp.gaps} ({round(100*hsp.gaps/alignment.length, 3)}%)")
                    print(f"  Number of positive (conservative) substitutions in the alignment: {hsp.positives} ({round(100*hsp.positives/alignment.length,3)}%)")
                    print(f"  Identities: {hsp.identities} ({round(100*hsp.identities/alignment.length,3)}%)\n")
                    print("Alignment:")
                    print(f"  Query:     {hsp.query}") 
                    print(f"  Match:     {hsp.match}")
                    print(f"  Subject:   {hsp.sbjct}\n")
                    if hsp.expect < 0.05:
                        print(f"Your E-value {hsp.expect} is lower than 0.05.\n")
                        print("An E-value lower than 0.05 is a good indicator of potential homology.\n")
                    elif hsp.expect >= 0.05:
                        print(f"Your E-value {hsp.expect} is equal or higher than 0.05\n")
                        print("You may consider your proteins as not homologues, please take into consideration.\n")
                    print("----------------------------------------")
    
    #save the same output as a text file:
    save_results_blastp(output_file, txt_name)

def save_results_blastp(output_file, txt_name):
    """
    
    This function takes the output_file from funtion run_blastp() and saves its results as txt file.
    """
    
    # parse the XML file using NCBIXML module
    with open(output_file, "r") as results:
        blast_records = NCBIXML.parse(results)
        
        # iterates over records (we only have one so we can use read - I need to ask you what you prefer)
        for record in blast_records:
            # iterates over alignments in the each record
            for alignment in record.alignments:
                #open file and save important information about results:
                f = open(txt_name, 'a')
                f.write(f"Alignment: {alignment.title}\n")
                f.write(f"Database searched: {record.database}\n")
                f.write(f"ID of protein searched: {record.query}\n")
                f.write(f"ID of protein matched: {alignment.hit_id}\n")
                if record.query.split(" ")[0] in protein_names:
                    f.write(f"ID {record.query.split(' ')[0]} is: '{protein_names[record.query.split(' ')[0]]}'\n")
                for hsp in alignment.hsps:
                    f.write(f"\n\n  E-value: {hsp.expect}\n")
                    f.write(f"  Bit Score: {hsp.bits}\n")
                    f.write(f"  Length: {alignment.length}\n")
                    f.write(f"  Number of gaps: {hsp.gaps} ({round(100*hsp.gaps/alignment.length,3)}%)\n")
                    f.write(f"  Identities: {hsp.identities} ({round(100*hsp.identities/alignment.length,3)}%)\n")
                    f.write(f"  Number of positive (conservative) substitutions in the alignment: {hsp.positives} ({round(100*hsp.positives/alignment.length,3)}%)\n")
                    f.write("\nAlignment:\n")
                    f.write(f"  Query:       {hsp.query}\n")
                    f.write(f"  Match:       {hsp.match}\n")
                    f.write(f"  Subject:     {hsp.sbjct}\n\n")
                    if hsp.expect < 0.05:
                        f.write(f"Your E-value {hsp.expect} is lower than 0.05,\n")
                        f.write("An E-value lower than 0.05 is a good indicator of potential homology.\n")
                    elif hsp.expect >= 0.05:
                        f.write(f"Your E-value {hsp.expect} is equal or higher than 0.05,\n")
                        f.write("You may consider your proteins as not homologues, please take into consideration.\n")
                    f.write("-------------------------------------------------------------\n")
                    f.close()
                print("\n")
        
################################################################################################################
# Option 2
################################################################################################################
def option_2():
    """
    This funtion prints insstruction for user for option 2 in main program.
    """
    print('\033[4m' + "\nPlease enter your fasta file name as followng:\n" + '\033[0m')
    print("\n~Please save your fasta file in the folder 'data' inside the program folder 'program_group6', as exlpained in manual.")
    print("~Make sure you file name is unique and ends with .faa or .fasta (for example: 'original_name.fasta').")
    print("~If a file with the same name already exists, it will delete it!\n")
    print("\033[95m" + "To go back to the main menu, enter 'end'" + "\033[0m")
    

##############################################################################################################################
# Option 3
##############################################################################################################################
    
def option_3():
    """
    This function is designed to facilitate the retrieval of a sequence based on its ID from a FASTA file.
    
    Returns:
    - dict: A dictionary containing sequences from the FASTA file.
    
    """

    # Call the functions for file path and create a dictionary:
    file_path = fasta_file_path()

    # Allow the user to go back to the main menu:
    if file_path == "end":
        print("\033[95m" + "\nYou chose to go back to the main menu.\n" + "\033[0m")
        return

    #Convert the FASTA file to a dictionary of sequences:
    sequences_dict = fasta_to_dict(file_path)

    #Ask the user for an ID, ensuring case insensitivity:
    while True:
        users_choice_id = input("\nWhat is the ID to look for in the file?\n" +
                                "\033[95m" + "To go back to the main menu, enter 'end'" + "\033[0m")

        #Check if the user-provided ID is in the dictionary:
        if users_choice_id in sequences_dict:
            print(f"ID: {users_choice_id}\nSequence: {sequences_dict[users_choice_id]}\n")
            break

        #Check if the provided ID, when case insensitive, matches any keys in the dictionary:
        for key in sequences_dict:
            if key.lower() == users_choice_id.lower():
                print(f"ID: {users_choice_id}\nSequence: {sequences_dict[key]}\n")
                break

        # let user option to go back to the main menu:
        if users_choice_id == "end":
            print("\033[95m" + "\nYou chose to go back to the main menu.\n" + "\033[0m")
            return
        else:
            print(f"Unfortunately, '{users_choice_id}' is not an ID. Please try again.\n")

            
            
def fasta_file_path():
    """
    This function is designed to open FASTA file and check if the name of the fasta was givven correctrly. 
    If not, it will let the user try again.
    
    Parameters:
    - file_path (str): The name to the FASTA file that is saved in "data" folder inside the program folder 'program_group6'.

    Returns:
    - dict: A dictionary where keys are sequence IDs and values are sequences.
    
    """
    while True:
        try:
            # Get the file path from user input. Input brings back a str.
            print('\033[4m' + "\nPlease enter you fasta file name as followng:\n" + '\033[0m')
            print("~Please save your fasta file in the folder 'data' inside the program file 'program_group6', ends with .faa or .fasta, as exlained in manual.")
            print("\033[95m" + "To go back to the main menu, enter 'end'" + "\033[0m")
            file_path = input("Enter the name of your FASTA file, do not use quotations: ")
            
            if file_path == "end":
                return file_path
            
            direct = "../data/" + file_path
            # Try to open the file to check if it exists
            with open(direct):
                #if the file opens, pass
                pass
            

            # If the file opening is successful, it will return the file path as str
            return direct
        
        except FileNotFoundError:
            print(f"Error: The file path '{file_path}' does not exist. Please try again.\n")
            

def fasta_to_dict(file_path):
        """
        This function takes the file path that we checked, and returns a dictionary of sequences.

        Parameters:
        - filepath (str): The path to the FASTA file.

        Returns:
        - dict: A dictionary where keys are sequence IDs and values are sequences.
        """
        sequences_dict = {}
        for record in SeqIO.parse(file_path, "fasta"):
            sequences_dict[record.id] = record.seq
        return sequences_dict


    
    
##############################################################################################################################
# Code for bash
##############################################################################################################################
#This part of the code runs our python functions that were called from the program (bash script).
#sys.argv is the arguments that were called in the bash script.
#sys.argv[1]=="adjust_csv" will call the function adjust_csv() with the csv saved as path "x".
#sys.argv[1]=="option_2" will call the function "option_2"
#sys.argv[1]=="option_3" will call the function "option_3"
#sys.argv[1] == "print_results_blastp" will call function print_results_blastp with arguments provided by user.
#sys.argv[1] == "Print_User" will run instructins for user (run from here for nicer visuallity)
#len(sys.argv) == 1 and sys.argv[0]== "program_functions.py" run the main menu - hence starting the program.

#To run our program: 
if __name__ == "__main__":
    if len(sys.argv) == 0:
        pass
    elif len(sys.argv) > 1:
        if sys.argv[1]=="adjust_csv":
            x = sys.argv[2]
            adjust_csv(x)
        if sys.argv[1]=="option_2":
            option_2()       
        elif sys.argv[1]=="option_3":
            option_3()
        elif sys.argv[1] == "print_results_blastp":
            print_results_blastp(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == "Print_User":
            print("Enter the Protein IDs you want to search for. If multiple protein IDs separate with OR (i.e. SOH98281.1 OR CAI77379.1).")
            print("Will only print results for valid IDs.")
            print( "\033[95m" + "To go back to the main menu, enter 'end'\n" + "\033[0m")
    elif len(sys.argv) == 1 and sys.argv[0]== "program_functions.py":
        user_main_menu()

    
