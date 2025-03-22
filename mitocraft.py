#!/usr/bin/env python

'''
Author - Devansh Shah
'''
import argparse
import os
import pandas as pd
import logging
import sys
import subprocess
import re
import random

# Custom exceptions for error handling
class MitocraftError(Exception):
    """Base exception for Mitocraft errors."""
    pass

class ReferenceBaseError(MitocraftError):
    """Exception raised when reference base doesn't match."""
    pass

class PipelineError(MitocraftError):
    """Exception raised when no suitable pipeline is found."""
    pass

class CommandError(MitocraftError):
    """Exception raised when a command execution fails."""
    pass

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='logging_main.log')
logger = logging.getLogger(__name__)

# MAGIC NUMBERS for the TALE-NT Tool
MIN_SPACER = 14
MAX_SPACER = 18
ARR_MIN = 14
ARR_MAX = 18
FILTER = 1 # keep 1 for specific efficient TALES ELSE keep 2 to generate all kinds of TALE! (see TALEN FAQ page)
CUT_POS = 31

def import_pipeline(pipeline_name):
    """Import the selected pipeline module from the pipelines folder."""
    try:
        module = __import__(f"pipelines.{pipeline_name}", fromlist=["process_mtDNA", "append_to_excel", "list_to_fasta"])
        process_mtDNA = module.process_mtDNA
        append_to_excel = module.append_to_excel
        list_to_fasta = module.list_to_fasta
    except ImportError as e:
        logger.error("Failed to import pipeline: %s", e)
        raise ValueError(f"Unknown pipeline: {pipeline_name}")

    return process_mtDNA, append_to_excel, list_to_fasta

def validate_input_files(input_file, additional_file):
    """Validate that the input files exist."""
    for file in [input_file, additional_file]:
        if not os.path.isfile(file):
            logger.error("The file - %s does not exist.", file)
            raise FileNotFoundError(f"The file - {file} does not exist.")

def run_findTAL(TALEN_dir, FASTA_fn, OUT_fn):
    """Run the findTAL script using subprocess."""
    findTAL_path = os.path.join(TALEN_dir, "findTAL.py")

    cmd = [
        "conda", "run", "-n", "run_talen_env",
        "python", findTAL_path,
        "--fasta", FASTA_fn,
        "--min", str(MIN_SPACER),
        "--max", str(MAX_SPACER),
        "--arraymin", str(ARR_MIN),
        "--arraymax", str(ARR_MAX),
        "--outpath", OUT_fn,
        "--filter", str(FILTER)
    ]

    # Add --filterbase only if FILTER is 1
    if FILTER == 1:
        cmd.append("--filterbase")
        cmd.append(str(CUT_POS))

    logger.info("[cmd] %s", ' '.join(cmd))  # Log the command for debugging
    try:
        result = subprocess.call(cmd)
        if result != 0:
            raise RuntimeError(f"Command failed with return code {result}.")
    except Exception as e:
        logger.error("Error executing command: %s", e)
        raise CommandError(f"Error executing command: {e}")


def process_talen_output(all_windows_df, bystanders_df, txt_file_path, output_excel_file_path):
    """Process the output from the TALEN tool and update the provided DataFrames."""
    logger.info("Loading TXT file from %s", txt_file_path)

    try:
        txt_df = pd.read_csv(txt_file_path, delimiter='\t', skiprows=2)
        logger.info("Successfully loaded TXT file.")
    except pd.errors.ParserError:
        logger.error("Error reading the TXT file. Ensure it is tab-delimited and correctly formatted.")
        raise
    except Exception as e:
        logger.error("Failed to load TXT file %s: %s", txt_file_path, e)
        raise

    if 'Plus strand sequence' not in txt_df.columns:
        logger.error("Column 'Plus strand sequence' not found in the TALEN file.")
        raise ValueError("Column 'Plus strand sequence' not found in the TALEN file")

    plus_strand_sequences = txt_df['Plus strand sequence'].tolist()
    txt_bases = set(re.findall(r'[a-z]+', ' '.join(plus_strand_sequences)))

    logger.info("Comparing Excel file sequences with TALEN_Tool TXT file bases.")
    for index, row in all_windows_df.iterrows():
        cleaned_sequence = row['Window Sequence'].translate(str.maketrans('', '', '{}[]'))
        all_windows_df.at[index, 'Matching TALEs'] = cleaned_sequence.lower() in txt_bases

    logger.info("Saving updated DataFrames to output Excel file %s.", output_excel_file_path)
    
    with pd.ExcelWriter(output_excel_file_path, engine='openpyxl') as writer:
        all_windows_df.to_excel(writer, sheet_name='All_Windows', index=False)
        bystanders_df.to_excel(writer, sheet_name='Bystander_Effects', index=False)

    logger.info("Updated final Excel file saved successfully.")

def append_match_info(all_windows_df, bystanders_df, txt_file_path, output_excel_file_path):
    logger.info("Loading TXT file for retrieving the TALE sequences from %s", txt_file_path)
    
    # Load the TXT file
    try:
        txt_df = pd.read_csv(txt_file_path, delimiter='\t', skiprows=2)
        logger.info("Successfully loaded TXT file.")
    except Exception as e:
        logger.error("Failed to load TXT file %s: %s", txt_file_path, e)
        raise

    if 'Plus strand sequence' not in txt_df.columns:
        logger.error("Column 'Plus strand sequence' not found in the TALEN file.")
        raise ValueError("Column 'Plus strand sequence' not found in the TALEN file")

    plus_strand_sequences = txt_df['Plus strand sequence'].tolist()
    txt_bases = set(re.findall(r'[a-z]+', ' '.join(plus_strand_sequences)))

    logger.info("Appending match information to the All_Windows DataFrame.")

    # Append results to the DataFrame
    for index, row in all_windows_df.iterrows():
        cleaned_sequence = row['Window Sequence'].translate(str.maketrans('', '', '{}[]')).lower()
        for txt_base in txt_bases:
            if txt_base == cleaned_sequence:  # Exact match check
                left_index = 1
                right_index = 1
                matching_sequences = [sequence for sequence in plus_strand_sequences if re.sub(r'[^a-z]', '', sequence) == txt_base]
                for sequence in matching_sequences:
                    lower_indices = [i for i, char in enumerate(sequence) if char.islower()]
                    if not lower_indices:
                        continue  # Skip if there are no lowercase letters
                    first_lower_index = lower_indices[0]
                    last_lower_index = lower_indices[-1]
                    L_part = sequence[:first_lower_index].strip()
                    R_part = sequence[last_lower_index + 1:].strip()
                    logger.info('index - %s, Left TALE - %s, Right TALE - %s', index, L_part, R_part)
                    all_windows_df.at[index, f'LeftTALE{left_index}'] = L_part
                    all_windows_df.at[index, f'RightTALE{right_index}'] = R_part
                    left_index += 1
                    right_index += 1

    logger.info("Done adding the TALE sequences!")

    # Only process bystanders if it's not empty
    if not bystanders_df.empty and 'Bystander Position' in bystanders_df.columns:
        logger.info("Filtering Bystander Effects to a maximum of three entries per Bystander Position.")
        bystanders_df = bystanders_df.groupby('Bystander Position').head(3).reset_index(drop=True)
        bystanders_df = bystanders_df.sort_values(by='Bystander Position').reset_index(drop=True)
    else:
        logger.info("No Bystander Data to process.")

    # Saving the updated DataFrames to the Excel file
    logger.info("Saving updated DataFrames to output Excel file %s.", output_excel_file_path)
    with pd.ExcelWriter(output_excel_file_path, engine='openpyxl') as writer:
        all_windows_df.to_excel(writer, sheet_name='All_Windows', index=False)
        if not bystanders_df.empty:  # Only write if not empty
            bystanders_df.to_excel(writer, sheet_name='Bystander_Effects', index=False)

    logger.info("Updated final Excel file saved successfully.")

def setup_directories(parent_directory):
    """Create necessary directories for output files."""
    running_directory = os.path.join(parent_directory, 'running')
    os.makedirs(running_directory, exist_ok=True)  # Create the 'running' directory

    directories = {
        'fasta': os.path.join(running_directory, 'fasta'),
        'pipeline_windows': os.path.join(running_directory, 'pipeline_windows'),
        'all_windows': os.path.join(running_directory, 'all_windows'),
        'talen': os.path.join(running_directory, 'talen'),
        'matching_output': os.path.join(running_directory, 'matching_output'),
        'final_output': os.path.join(parent_directory, 'final_output')  # Changed to parent_directory
    }

    for dir_path in directories.values():
        os.makedirs(dir_path, exist_ok=True)
        
    return directories

def fasta_to_txt(fasta_file, txt_file):
    """Convert a FASTA file to a plain text file containing only the sequence."""
    with open(fasta_file, "r") as fh:
        # Join lines that do not start with '>' into a single sequence
        sequence = ''.join(line.strip() for line in fh if not line.startswith('>'))
    
    # Write only the sequence to the TXT file
    with open(txt_file, "w") as out_fh:
        out_fh.write(sequence)
    logger.info(f"Converted {fasta_file} to {txt_file} successfully.")


def main():
    parser = argparse.ArgumentParser(description='Process DNA sequence for base editing.')
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Set the default input file path
    default_input_file = os.path.join(script_dir, 'inputs', 'mito.txt')

    parser.add_argument('--input_file', type=str, default=default_input_file, help='File containing the mtDNA sequence')
    #parser.add_argument('--input_file', type=str, default='inputs/mito.txt', help='File containing the mtDNA sequence') #hard code?
    parser.add_argument('position', type=int, help='Position of the base to be changed')
    parser.add_argument('reference_base', type=str, help='Original base of the position')
    parser.add_argument('mutant_base', type=str, help='Mutant base to be changed into')
    args = parser.parse_args()

    # Convert reference and mutant bases to uppercase
    reference_base = args.reference_base.upper()
    mutant_base = args.mutant_base.upper()

    # Define file paths
    input_file = args.input_file
    #additional_file = 'inputs/annotated_human_mtDNA_10022024_for_bystanders.xlsx' #THIS IS ONLY FOR THE HUMAN MITOCHONDRIAL GENOME!!
    # additional_file = os.path.join(script_dir, 'inputs', 'annotated_human_mtDNA_10022024_for_bystanders.xlsx')
    additional_file = os.path.join(script_dir, 'inputs', 'annotated_human_mtDNA_10022024_for_bystanders_EDITED.xlsx') # this ONLY includes the mutations possible by the base editors

    # Validate input files
    validate_input_files(input_file, additional_file)

    # Read the input DNA sequence
    logger.info("Reading the input DNA sequence %s.", input_file)

    # Check if the input file is a FASTA file by its extension
    if input_file.endswith('.fasta') or input_file.endswith('.fa'):
        txt_file = input_file.rsplit('.', 1)[0] + '.txt'  # Create a corresponding .txt file name
        fasta_to_txt(input_file, txt_file)  # Convert FASTA to TXT
        with open(txt_file, "r") as fh:
            mtDNA_seq = fh.read().replace("\n", "")
            logger.info("main txt file %s", mtDNA_seq)
    else:
        with open(input_file, "r") as fh:
            mtDNA_seq = fh.read().replace("\n", "")

    # Check if the reference base matches the input file base at the specified position
    if reference_base != mtDNA_seq[args.position - 1].upper():
        logger.error("Incorrect reference base at position %d. User Input: %s, Actual Reference: %s", 
                     args.position, reference_base, mtDNA_seq[args.position - 1].upper())
        print(f"Incorrect reference base at position {args.position}. Expected: {reference_base}, Found: {mtDNA_seq[args.position - 1].upper()}")
        raise ReferenceBaseError(f"Incorrect reference base at position {args.position}. Expected: {reference_base}, Found: {mtDNA_seq[args.position - 1].upper()}")

     # Define pipelines based on reference and mutant bases
    if (reference_base, mutant_base) in [('C', 'T'), ('G', 'A')]:
        pipelines = ["Mok2020_G1397", "Mok2020_G1333", "Mok2022_DddA11"]
        print("Not editable by the Cho pipeline")
    elif (reference_base, mutant_base) in [('A', 'G'), ('T', 'C')]:
        pipelines = ["Cho_sTALEDs"]
        print("Not editable by the Mok pipelines")
    else:
        print(f"No pipeline found for reference base {reference_base} and the mutant base {mutant_base}")
        logger.error("No pipeline found for reference base %s and mutant base %s.", reference_base, mutant_base)
        raise PipelineError(f"No pipeline found for reference base {reference_base} and the mutant base {mutant_base}")

    all_windows_combined = pd.DataFrame()  
    all_bystanders_combined = pd.DataFrame()

    # Process each pipeline
    for pipeline in pipelines:
        logger.info("Processing pipeline: %s", pipeline)

        # Import the selected pipeline functions
        process_mtDNA, append_to_excel, list_to_fasta = import_pipeline(pipeline)

        logging.info("Processing mtDNA sequence for position %d using pipeline %s.", args.position, pipeline)
        all_windows, adjacent_bases = process_mtDNA(mtDNA_seq, args.position)

        if not adjacent_bases:
            logger.warning("The base found at position %d cannot edited by the pipeline %s.", args.position, pipeline)
            continue

        # Define paths for output files
        parent_directory = os.path.dirname(os.path.abspath(__file__))
        directories = setup_directories(parent_directory)

        # Define output file paths
        fasta_file = os.path.join(directories['fasta'], f'adjacent_bases_{args.position}.fasta')
        allw_file = os.path.join(directories['pipeline_windows'], f'{pipeline}_{args.position}.xlsx')

        # Write adjacent bases to FASTA file
        logging.info("Writing adjacent bases to FASTA file.")
        fasta_content = list_to_fasta(adjacent_bases, args.position)
        with open(fasta_file, 'w') as file:
            file.write(fasta_content)
            logging.info("Finished writing FASTA file to %s.", fasta_file)

        # Append all windows to the combined DataFrame
        logger.info("Writing all windows to Excel file.")
        windows_df, bystanders_df = append_to_excel(all_windows, additional_file, allw_file)

        # Concatenate to combined DataFrames
        all_windows_combined = pd.concat([all_windows_combined, windows_df], ignore_index=True)
        all_bystanders_combined = pd.concat([all_bystanders_combined, bystanders_df], ignore_index=True)

     # After processing all pipelines, save both combined DataFrames to one Excel file
    output_combined_path = os.path.join(directories['all_windows'], f'all_windows_{args.position}.xlsx')
    logger.info("Saving combined outputs to %s.", output_combined_path)

    with pd.ExcelWriter(output_combined_path, engine='openpyxl') as writer:
        all_windows_combined.to_excel(writer, sheet_name='All_Windows', index=False)
        all_bystanders_combined.to_excel(writer, sheet_name='Bystander_Effects', index=False)

    logger.info("Combined output saved successfully.")

    # Now call the TALEN tool
    #TALEN_dir = "software/talent_tools_master"
    TALEN_dir = os.path.join(script_dir, "software", "talent_tools_master")
    out_fn = os.path.join(directories['talen'], f'TALENT_{args.position}.txt')
    run_findTAL(TALEN_dir, fasta_file, out_fn)
    
    # Process the TALEN output
    output_excel_path = os.path.join(directories['matching_output'], f'matching_tales_{args.position}.xlsx')
    process_talen_output(all_windows_combined, all_bystanders_combined , out_fn, output_excel_path)

    f_out = os.path.join(directories['final_output'], f'final_{args.position}.xlsx')
    #append_match_info(all_windows_combined, all_bystanders_combined, out_fn, f_out)
     # Check if user uploaded a file
    #if args.input_file == 'inputs/mito.txt':  # User did not upload an input file
    if args.input_file == os.path.join(script_dir, 'inputs', 'mito.txt'):
        append_match_info(all_windows_combined, all_bystanders_combined, out_fn, f_out)
    else:
        all_bystanders_combined = pd.DataFrame()  # Clear bystanders DataFrame
        append_match_info(all_windows_combined, all_bystanders_combined, out_fn, f_out)

if __name__ == "__main__":
    main()
