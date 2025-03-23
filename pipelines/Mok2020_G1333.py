#!/usr/bin/env python

'''
Author - Devansh Shah
'''
PIPELINE = "Mok2020_G1333"

#importing packages
import argparse
import random 
import os 
import pandas as pd
import logging

#setting up the logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='logging_Mok2020_G1333.log')
logger = logging.getLogger(__name__)

#MARKING THE DNA SEQUENCE FOR THE TARGET BASE AND BYSTANDER
# Target base with square brackets [] : Off-target base with curly braces {}
def mark_bases(sequence, target_position, off_target_positions):
    """Mark the target and bystander bases in the window"""
    logger.debug("Marking bases in the sequence.")
    target_position -= 1 #one indexed
    off_target_positions = set(p - 1 for p in off_target_positions)
    marked_sequence = [] #making an empty list
    for index, char in enumerate(sequence):
        if index == target_position:
            marked_sequence.append(f"[{char}]")  # Target base with square brackets []
        elif index in off_target_positions:
            marked_sequence.append(f"{{{char}}}")  # Off-target base with curly braces {}
        else:
            marked_sequence.append(char)  # No special marking
    return ''.join(marked_sequence)

def mark_base_at_position(sequence, target_position):
    """Mark the base at the target position --> to mark the ADJACENT off-targets"""
    logger.debug("Marking base at the specific position")
    marked_sequence = []
    for index, char in enumerate(sequence):
        if index == target_position:
            marked_sequence.append(f"{{{char}}}")  # Target ADJACENT base with curly brackets {}
        else:
            marked_sequence.append(char)  # No special marking
    return ''.join(marked_sequence)

def reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence"""
    logger.debug("Generating reverse complement.")
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': random.choice(['A', 'T', 'C', 'G']), '}' : '{', '{' : '}', '[':']', ']':'['} #why have I added the complement of the brackets too here?
    reverse_complement_sequence = ''.join([complement_dict.get(base, base) for base in sequence])
    return reverse_complement_sequence[::-1]

def complementing(sequence):
    """Get the complement of a DNA sequence."""
    logger.debug("Generating complement.")
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': random.choice(['A', 'T', 'C', 'G'])}
    complement_sequence = ''.join([complement_dict.get(base, base) for base in sequence])
    return complement_sequence

def create_window(mtDNA_seq, pos, start_index, end_index):
    """To create the window"""
    logger.debug("Creating window from position %d.", pos)
    window = mtDNA_seq[start_index:end_index]
    return window

def generate_TC_windows(mtDNA_seq, pos, window_size):
    """Generating windows in a 5'TC-context from the 5'end."""
    logger.debug("Generating TC windows for position %d with window size %d from the 5'-end.", pos, window_size)
    TC_windows = [] 
    for i in range(4, 12): 
        if i != 11: 
            start_index = max(0, pos - i)  # starting index of window from the 5'end (5'-3' TC context)
            end_index = min(len(mtDNA_seq), start_index + window_size)  # ending index for the window
            if len(mtDNA_seq[start_index:end_index]) == window_size: #to check whether the window size is valid!
                window = create_window(mtDNA_seq, pos, start_index, end_index) #creating the window
                TC_windows.append(window) #adding the generated window to the list
    return TC_windows

def generate_GA_windows(mtDNA_seq, pos, window_size):
    """Generating windows in a 5'GA-context from the 3' end."""
    logger.debug("Generating GA windows for position %d with window size %d from the 3' end.", pos, window_size)
    GA_windows = [] #storing the generated windows in this list
    for i in range(4, 12):  # loop over position 4-10 from the left
        if i != 11: #excluding the 11th pos! 
            start_index = max(0, pos - window_size + i -1) #start index of window from the 3' end (since we are looking at 3'-5' GA context)
            end_index = min(len(mtDNA_seq), start_index + window_size)  # end index for the window
            if len(mtDNA_seq[start_index:end_index]) == window_size: #checking whether the window is valid or not
                window = create_window(mtDNA_seq, pos, start_index, end_index)
                GA_windows.append(window[::-1]) # reversing to read from 5'-3' direction!! (imp in terms of reading sequences!)
    return GA_windows

def remove_whitespace(seq):
    """Remove any space in input sequence"""
    return ''.join(seq.split())

def capitalize(seq):
    """Capitalize the input sequence"""
    return seq.upper()

def lowerize(seq):
    """Lowerize the input sequence"""
    return seq.lower()

def find_consecutive_GA_sequences(sequence):
    """Find 5'-GA contexts in the entire sequence"""
    logger.debug("Finding 5'-GA contexts.")
    GA_positions = []
    start_index = 0
    while True:
        index = sequence.find('GA', start_index)
        if index == -1:
            break
        GA_positions.append(index + 1)
        start_index = index + 1
    return GA_positions

def find_consecutive_TC_sequences(sequence):
    """Find 5'-TC contexts in the entire sequence"""
    logger.debug("Finding 5'-TC sequences.")
    TC_positions = []
    start_index = 0
    while True:
        index = sequence.find('TC', start_index)
        if index == -1:
            break
        TC_positions.append(index + 1 + 1) # +1 because i want to store the position of the C in 5'-TC
        start_index = index + 1
    return TC_positions

def count_GA_sequences(sequence):
    """Find all 5'-GA contexts in the window"""
    logger.debug("Counting all the 5'-GA contexts in the window.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('GA', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_TC_sequences(sequence):
    """Find all 5'-TC contexts in the window"""
    logger.debug("Counting all the 5'-TC contexts in the window.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('TC', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def list_to_fasta(dna_list, pos):
    """Converting the adjacent_bases into FASTA format"""
    logger.debug("converting to FASTA format")
    fasta_str = ""
    sequence = dna_list
    header = f">chrM_{pos}"
    fasta_str = f"{header}\n{sequence}\n"
    return fasta_str

def find_GA_positions(window, start_position):
    """Find positions of 5'-GA contexts in the window relative to the mtDNA sequence."""
    logger.debug("Finding 5'-GA contexts in the window.")
    ga_positions = []
    start_index = 0
    while True:
        index_ga = window.find('GA', start_index)
        if index_ga == -1:
            break
        ga_positions.append(start_position + index_ga + 1)  # +1 for 1-indexed
        start_index = index_ga + 1
    return ga_positions

def find_TC_positions(window, start_position):
    """Find positions of 5'-TC contexts in the window relative to the mtDNA sequence."""
    logger.debug("Finding 5'-TC contexts in the window.")
    tc_positions = [] 
    start_index = 0
    while True:
        index_tc = window.find('TC', start_index)
        if index_tc == -1:
            break
        tc_positions.append(start_position + index_tc + 1 + 1)  # +1 for 1-indexed and another +1 to include position of the C
        start_index = index_tc + 1
    return tc_positions

def process_mtDNA(mtDNA_seq, pos):
    """Main function which processes the DNA sequence."""
    logger.info("Processing DNA sequence for position %d.", pos)

    #processing the DNA sequence 
    nospace_mtDNA = capitalize(remove_whitespace(mtDNA_seq))
    consecutive_TC_positions = find_consecutive_TC_sequences(nospace_mtDNA) 
    consecutive_GA_positions = find_consecutive_GA_sequences(nospace_mtDNA)
    dummy = 0
    
    if pos in consecutive_TC_positions: #found a 5'-TC context
        logger.info("Base at position %d is in a 5'-TC context.", pos)
        ref, mut, all_windows, dum = 'C', 'T', [], []

        circular_seq = nospace_mtDNA + nospace_mtDNA # circularizing the DNA sequence (to keep the mtDNA sequence functional - good thing that the D-Loop is present)
        
        start_index = pos - (16 + 15) #(target - 30) 
        end_index = pos + (15 + 15) # (target + 30)

        adjacent_bases = circular_seq[start_index:end_index] #splicing the base sequence around the target
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases))
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)

        FLAG=None #this is to check if the EXACT right adjacent base is also a C or not. how should I look if SECOND adjacent base also is part of the same codon??
        if right_adjacent_bases[0] == 'C':
            dummy = 1 # to start counting for bystander edits
            dum.append(pos+1) #DNA pos of the right adjacent base
            FLAG=True

        for window_size in range(14, 19): #ws=14bp, 15bp, 16bp, 17bp, 18bp
            TC_windows = generate_TC_windows(circular_seq, pos, window_size) #since it is in a 5'-TC context, will generate windows from the 5'-end!
            TALES = False #this is used later for (TALE-NT tool)

            for num, window in enumerate(TC_windows, start=4):
                window_description = f"Position {num} from the 5' end" #from 5' end because we are generating windows from that end
                ws = f"{window_size}bp"

                off_target_sites = count_TC_sequences(window[2:10]) + count_GA_sequences(window[-10:-2]) + dummy #count of ALL BYSTANDER EDITS

                #marking the above BYSTANDERS (5'-TC from the 5' end or 5'-GA from the 3'end) at positions 4-10 from either end
                marked_window = mark_bases(window, num, [(x + window_size - 10) for x in find_consecutive_GA_sequences(window[-10:-2])] + [(x + 2) for x in find_consecutive_TC_sequences(window[2:10])])
                
                # Store the positions of `{` and `}` in separate lists
                brace_left_positions = [i for i, char in enumerate(marked_window) if char == '{']
                brace_right_positions = [i for i, char in enumerate(marked_window) if char == '}']

                ftc, fga = [],[]
                
                #finding the position of the bystanders (5'-TC context from the 5'-end) and (5'-GA context from the 3'-end)
                start_position = pos - num + 1 #this adjusts the start_position to the first position from the 5'-end
                ga_positions = find_GA_positions(window[-10:-2], start_position + window_size - 10 - 1)
                tc_positions = find_TC_positions(window[2:10], start_position + 2 - 1)
                ftc = tc_positions + dum #adding the position of the (5'-TCC) context too NEXT TO TARGET!
                fga = ga_positions
                
                if pos in ftc: #removing the position of the target itself!
                    ftc.remove(pos)
                    #off_target_sites=+1
                
                for left_pos, right_pos in zip(brace_left_positions, brace_right_positions):
                    # Initial offset
                    offset_tc, offset_ga = 0, 0

                    # Adjust the offset based on the position in brace_left_positions
                    if left_pos == brace_left_positions[0]:  # First item in the brace_left_positions
                        offset_tc = 2  # For 5'-T{C}C context, offset by 2 for the first match
                    else:  # For subsequent items
                        offset_tc += 2  # Increment offset for subsequent matches

                        # Adjust the offset based on the position in brace_left_positions
                    if left_pos == brace_left_positions[0]:  # First item in the brace_left_positions
                        offset_ga = 0  # For 5'-{G}GA context, offset by 2 for the first match
                    else:  # For subsequent items
                        offset_ga += 2  # Increment offset for subsequent matches

                    # Check if the left_pos + 1 < num (the condition for the [ ] brackets)
                    if left_pos + 1 > num:
                        offset_tc += 2  # If true, we add an offset of 2
                        offset_ga += 2

                    # Check for 5'-TC{C} context: 'T' before '{' and 'C' after '}'
                    if left_pos > 0 and marked_window[left_pos - 1] == 'T' and marked_window[right_pos + 1] == 'C':
                        # Mark the second C (right_pos + 1)
                        mark_pos = right_pos + 1
                        marked_window = mark_base_at_position(marked_window, mark_pos)
                        off_target_sites += 1  # Increment off-target sites
                        ftc.append(mark_pos-offset_tc+start_position)
                        
                    # Check for 5'-{G}GA context: 'G' before '{' and 'A' after '}'
                    if right_pos < len(marked_window) - 1 and marked_window[left_pos - 1] == 'G' and marked_window[right_pos + 1] == 'A':
                        # Mark the first G (left_pos - 1)
                        mark_pos = left_pos - 1
                        marked_window = mark_base_at_position(marked_window, mark_pos)
                        off_target_sites += 1  # Increment off-target sites
                        fga.append(mark_pos+start_position-offset_ga)  # Store the position of first G

                 # Now handle marking for 5'-T[C]{C} context outside the loop
                int_pos = marked_window.find(']')
                if int_pos != -1 and marked_window[int_pos + 1] == 'C':  # Confirm it's a 5'-T[C]{C}
                    mark_pos = int_pos + 1
                    marked_window = mark_base_at_position(marked_window, mark_pos)
                        
                # Convert the final marked window to a string
                final_window_str = marked_window
                
                sorted_positions = sorted(ftc + fga)   #sorted version of the positions
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window_str, window_description, off_target_sites-1, sorted_positions, TALES, FLAG))
    
    elif pos in consecutive_GA_positions: #found a 5'-GA context
        logger.info("Base at position %d is in a 5'-GA context.", pos)
        ref, mut, all_windows, dum = 'G', 'A', [], []

        comple_nospace_mtDNA = complementing(nospace_mtDNA)
        circular_seq = nospace_mtDNA + nospace_mtDNA

        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)

        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases))
        
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", reverse_complement(left_adjacent_bases[::-1]), reverse_complement(right_adjacent_bases[::-1]))
        FLAG=None

        if left_adjacent_bases[-1] == 'G':
            dummy = 1
            dum.append(pos-1)
            FLAG=True

        for window_size in range(14, 19):
            GA_windows = generate_GA_windows(comple_nospace_mtDNA + comple_nospace_mtDNA, pos, window_size)
            TALES = False

            for num, window in enumerate(GA_windows, start=4):
                window_description = f"Position {num} from the 3' end"
                ws = f"{window_size}bp"

                off_target_sites = count_TC_sequences(window[2:10]) + count_GA_sequences(window[-10:-2]) + dummy

                reverse_window = complementing(window[::-1]) #to get the sequence in the top strand
                marked_window = mark_bases(reverse_window, window_size - num + 1,[(x + window_size - 10) for x in find_consecutive_GA_sequences(reverse_window[-10:-2])] + [(x + 2) for x in find_consecutive_TC_sequences(reverse_window[2:10])])
                
                # Store the positions of `{` and `}` in separate lists
                brace_left_positions = [i for i, char in enumerate(marked_window) if char == '{']
                brace_right_positions = [i for i, char in enumerate(marked_window) if char == '}']

                ftc, fga = [],[]
                
                start_position = pos - (window_size - num)
                tc_positions = find_TC_positions(reverse_window[2:10], start_position + 2 - 1)
                ga_positions = find_GA_positions(reverse_window[-10:-2], start_position +  window_size - 10 - 1) 
                ftc = tc_positions
                fga = ga_positions + dum

                if pos in fga:
                    fga.remove(pos)

                for left_pos, right_pos in zip(brace_left_positions, brace_right_positions):
                    # Initial offset
                    offset_tc, offset_ga = 0, 0

                    # Adjust the offset based on the position in brace_left_positions
                    if left_pos == brace_left_positions[0]:  # First item in the brace_left_positions
                        offset_tc = 2  # For 5'-T{C}C context, offset by 2 for the first match
                    else:  # For subsequent items
                        offset_tc += 2  # Increment offset for subsequent matches

                        # Adjust the offset based on the position in brace_left_positions
                    if left_pos == brace_left_positions[0]:  # First item in the brace_left_positions
                        offset_ga = 0  # For 5'-{G}GA context, offset by 0 for the first match
                    else:  # For subsequent items
                        offset_ga += 2  # Increment offset for subsequent matches

                    # Check if the left_pos + 1 > num (the condition for the [ ] brackets)
                    if left_pos + 1 > num:
                        offset_tc += 2  # If true, we add an offset of 2
                        offset_ga += 2

                    # Check for 5'-TC{C} context: 'T' before '{' and 'C' after '}'
                    if left_pos > 0 and marked_window[left_pos - 1] == 'T' and marked_window[right_pos + 1] == 'C':
                        # Mark the second C (right_pos + 1)
                        mark_pos = right_pos + 1
                        marked_window = mark_base_at_position(marked_window, mark_pos)
                        off_target_sites += 1  # Increment off-target sites
                        #ftc.append(mark_pos+start_position-offset_tc)  # Store the position of second C
                        ftc.append(mark_pos-offset_tc+start_position)
                        
                    # Check for 5'-{G}GA context: 'G' before '{' and 'A' after '}'
                    if right_pos < len(marked_window) - 1 and marked_window[left_pos - 1] == 'G' and marked_window[right_pos + 1] == 'A':
                        # Mark the first G (left_pos - 1)
                        mark_pos = left_pos - 1
                        marked_window = mark_base_at_position(marked_window, mark_pos)
                        off_target_sites += 1  # Increment off-target sites
                        fga.append(mark_pos+start_position-offset_ga)  # Store the position of first G
                
                int_pos = marked_window.find('[')
                if int_pos != -1 and marked_window[int_pos - 1] == 'G':  # Confirm it's a 5'-{G}[G]A
                    mark_pos = int_pos - 1
                    marked_window = mark_base_at_position(marked_window, mark_pos)
                # Convert the final marked window to a string
                final_window_str = marked_window
                
                sorted_positions = sorted(ftc + fga)  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window_str, window_description, off_target_sites-1, sorted_positions, TALES, FLAG))
    
    else:
        logger.warning("Base at position %d is not in a editable context and cannot be edited by the %s pipeline.", pos, PIPELINE)
        print(f"Position {pos} is not in a editable context and cannot be edited by the {PIPELINE}.")
        #logger.info("%s", nospace_mtDNA)
        return [], []  # Return empty lists to indicate failure

    return all_windows, adjacent_bases

def append_to_excel(all_windows, additional_file, output_file):
    """Search positions from 'ftc+fga' in the annotated Excel file and append to all_windows."""
    logger.info("Appending additional bystanders information to the Excel file.")
    
    # Create a DataFrame from all_windows
    all_windows_df = pd.DataFrame(all_windows, columns=[
        'Pipeline', 'Position', 'Reference Base', 'Mutant Base', 'Window Size', 
        'Window Sequence', 'Target Location', 'Number of Bystanders', 
        'Position of Bystanders', 'Optimal Flanking TALEs', 'Flag-CheckBystanderEffect'
    ])

    # Only read and process the additional file if it is provided
    if additional_file:
        if not os.path.isfile(additional_file):
            logger.warning("The additional bystander file - %s does not exist. Skipping appending bystander information.", additional_file)
            new_data = pd.DataFrame()  # Create an empty DataFrame
        else:
            ftc_fga_positions = set(pos for _, _, _, _, _, _, _, _, positions, _, _ in all_windows for pos in positions)
            additional_df = pd.read_excel(additional_file)
            filtered_df = additional_df[additional_df['mtDNA_pos'].isin(ftc_fga_positions)]
            new_data = filtered_df[['mtDNA_pos', 'Ref. Allele', 'Mutant Allele', 
                                     'Location', 'Predicted Impact', 'Syn vs NonSyn', 
                                     'AA Variant', 'Func. Impact', 'MutationAssessor Score']]
            new_data.columns = ['Bystander Position', 'Reference Base', 'Mutant Base', 
                                'Location On Genome', 'Predicted Mutation Impact', 
                                'SNV Type', 'AA Variant', 'Functional Impact', 
                                'MutationAssessor Score']
    else:
        logger.warning("No additional file provided. Skipping bystander information.")
        new_data = pd.DataFrame()  # Create an empty DataFrame if no additional file is provided

    with pd.ExcelWriter(output_file, engine='openpyxl', mode='a' if os.path.exists(output_file) else 'w') as writer:
        # Write all_windows to a separate sheet
        all_windows_df.to_excel(writer, sheet_name='All_Windows', index=False)
        # Append new data to a separate sheet only if it has data
        if not new_data.empty:
            new_data.to_excel(writer, sheet_name='Bystanders_Info', index=False)

    logger.info("Successfully appended bystander information to the Excel file, if available.")
    
    # Return the all_windows DataFrame for concatenation
    return all_windows_df, new_data


def main():
    parser = argparse.ArgumentParser(description='Process mtDNA sequence for base editing.')
    parser.add_argument('input_file', type=str, help='File containing the mtDNA sequence')
    parser.add_argument('position', type=int, help='Position of the base to be changed (between 1 and 16569)')
    parser.add_argument('additional_file', type=str, help='Excel file containing additional bystander information')
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        logger.error("The mtDNA file - %s does not exist.", args.input_file)
        return
    
    if not os.path.isfile(args.additional_file):
        logger.error("The additional bystander file - %s does not exist.", args.additional_file)
        return
    
    logger.info("Reading the input DNA sequence %s.", args.input_file)
    with open(args.input_file, "r") as fh:
        mtDNA_seq = fh.read().replace("\n", "")

    while True:  # Added retry loop
        logger.info("Processing mtDNA sequence for position %d.", args.position)
        all_windows, adjacent_bases = process_mtDNA(mtDNA_seq, args.position)

        # Check if editing is possible
        if not adjacent_bases:
            retry = input("Would you like to try a different position? (y/n): ").strip().lower()
            if retry == 'y':
                new_position = input("Enter a new position (between 1 and 16569): ")
                try:
                    new_position = int(new_position)
                    if new_position < 1 or new_position > 16569:
                        print("Position must be between 1 and 16569.")
                        continue  # Prompt for a new position
                    args.position = new_position  # Update the position
                    continue  # Restart the loop with the new position
                except ValueError:
                    print("Invalid input. Please enter a valid integer.")
                    continue  # Prompt for a new position
            else:
                logger.info("Exiting the program.")
                return  # Exit the program

        # Define paths for output files
        parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        fasta_directory = os.path.join(parent_directory, f'{PIPELINE}_fasta')
        csv_directory = os.path.join(parent_directory, f'{PIPELINE}_all_windows')

        # Defining the files
        file_name = f'{PIPELINE}_adjacent_bases_{args.position}.fasta'
        allw_name = f'{PIPELINE}_all_windows_{args.position}.xlsx'

        file_path = os.path.join(fasta_directory, file_name)
        allw_path = os.path.join(csv_directory, allw_name)

        # Making directory if it doesn't exist
        os.makedirs(fasta_directory, exist_ok=True)
        os.makedirs(csv_directory, exist_ok=True)

        if adjacent_bases:
            logger.info("Writing adjacent bases to FASTA file.")
            fasta_content = list_to_fasta(adjacent_bases, args.position)
            with open(file_path, 'w') as file:
                file.write(fasta_content)
                logger.info("Finished writing FASTA file to %s.", file_path)
        else:
            logger.info(f"Not possible to edit this base using the current pipeline")

        if all_windows:
            logger.info("Writing all windows and the bystander information to Excel file.")
            append_to_excel(all_windows, args.additional_file, allw_path)
        
        break  # Exit the retry loop if processing was successful

if __name__ == "__main__":
    main()
