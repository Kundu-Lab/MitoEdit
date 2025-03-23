#!/usr/bin/env python

'''
Author - Devansh Shah
'''
PIPELINE = "Mok2022_G1397_DddA11"
import argparse
import random 
import os 
import pandas as pd
import logging

#setting up the logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='logging_Mok2022_G1397_DddA11_with_Bystanders.log')
logger = logging.getLogger(__name__)

# MARKING THE DNA SEQUENCE FOR THE TARGET BASE AND BYSTANDER
def mark_bases(sequence, target_position, off_target_positions):
    """Mark the target and bystander bases in the window"""
    logger.debug("Marking bases in the sequence.")
    target_position -= 1
    off_target_positions = set(p - 1 for p in off_target_positions)
    marked_sequence = []
    for index, char in enumerate(sequence):
        if index == target_position:
            marked_sequence.append(f"[{char}]")  # Target base with square brackets []
        elif index in off_target_positions:
            marked_sequence.append(f"{{{char}}}")  # Off-target base with curly braces {}
        else:
            marked_sequence.append(char)  # No special marking
    return ''.join(marked_sequence)

def mark_base_at_position(sequence, target_position):
    """Mark the base at the target position --> to mark the target base in the window"""
    logger.debug("Marking base at the specific position")
    marked_sequence = []
    for index, char in enumerate(sequence):
        if index == target_position:
            marked_sequence.append(f"{{{char}}}")  # Target base with curly brackets {}
        else:
            marked_sequence.append(char)  # No special marking
    return ''.join(marked_sequence)

def reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence"""
    logger.debug("Generating reverse complement.")
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': random.choice(['A', 'T', 'C', 'G']), '}' : '{', '{' : '}', '[':']', ']':'['}
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

def generate_NC_windows(mtDNA_seq, pos, window_size): #THIS WILL BE --> AC,CC,TC
    """Generating windows in a NC-conext"""
    logger.debug("Generating NC windows for position %d with window size %d.", pos, window_size)
    NC_windows = []
    for i in range(4, 9):
        if i != 8:
            start_index = max(0, pos - window_size + i - 1)
            end_index = min(len(mtDNA_seq), start_index + window_size)
            if len(mtDNA_seq[start_index:end_index]) == window_size:
                window = create_window(mtDNA_seq, pos, start_index, end_index)
                NC_windows.append(window)
    return NC_windows

def generate_GN_windows(mtDNA_seq, pos, window_size): #THIS WILL BE --> GA, GG, GT
    """Generating windows in a GN-conext"""
    logger.debug("Generating GN windows for position %d with window size %d.", pos, window_size)
    GN_windows = []
    for i in range(4, 9):
        if i != 8:
            start_index = max(0, pos - i)
            end_index = min(len(mtDNA_seq), start_index + window_size)
            if len(mtDNA_seq[start_index:end_index]) == window_size:
                window = create_window(mtDNA_seq, pos, start_index, end_index)
                GN_windows.append(window[::-1]) #to append in the correct direction (but still on the bottom strand)
    return GN_windows

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
    """Find GA contexts in the entire sequence"""
    logger.debug("Finding consecutive GA sequences.")
    GA_positions = []
    start_index = 0
    while True:
        index = sequence.find('GA', start_index)
        if index == -1:
            break
        GA_positions.append(index + 1)
        start_index = index + 1
    return GA_positions

def find_consecutive_GT_sequences(sequence):
    """Find GT contexts in the entire sequence."""
    logger.debug("Finding consecutive GT sequences.")
    GT_positions = []
    start_index = 0
    while True:
        index = sequence.find('GT', start_index)
        if index == -1:
            break
        GT_positions.append(index + 1)  # 1-indexed position of the G base
        start_index = index + 1
    return GT_positions

def find_consecutive_GG_sequences(sequence):
    """Find GG contexts in the entire sequence."""
    logger.debug("Finding consecutive GG sequences.")
    GG_positions = []
    start_index = 0
    while True:
        index = sequence.find('GG', start_index)
        if index == -1:
            break
        GG_positions.append(index + 1)  # 1-indexed position of the first G
        start_index = index + 1
    return GG_positions


def find_consecutive_TC_sequences(sequence):
    """Find TC contexts in the entire sequence"""
    logger.debug("Finding consecutive TC sequences.")
    TC_positions = []
    start_index = 0
    while True:
        index = sequence.find('TC', start_index)
        if index == -1:
            break
        TC_positions.append(index + 1 + 1)
        start_index = index + 1
    return TC_positions

def find_consecutive_AC_sequences(sequence):
    """Find AC contexts in the entire sequence."""
    logger.debug("Finding consecutive AC sequences.")
    AC_positions = []
    start_index = 0
    while True:
        index = sequence.find('AC', start_index)
        if index == -1:
            break
        AC_positions.append(index + 2)  # 1-indexed position of the C base
        start_index = index + 1
    return AC_positions

def find_consecutive_CC_sequences(sequence):
    """Find CC contexts in the entire sequence."""
    logger.debug("Finding consecutive CC sequences.")
    CC_positions = []
    start_index = 0
    while True:
        index = sequence.find('CC', start_index)
        if index == -1:
            break
        CC_positions.append(index + 2)  # 1-indexed position of the second C base
        start_index = index + 1
    return CC_positions

def count_GA_sequences(sequence):
    """Count occurrences of GA in the sequence."""
    logger.debug("Counting GA sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('GA', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_GT_sequences(sequence):
    """Count occurrences of GT in the sequence."""
    logger.debug("Counting GT sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('GT', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_GG_sequences(sequence):
    """Count occurrences of GG in the sequence."""
    logger.debug("Counting GG sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('GG', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_TC_sequences(sequence):
    """Count occurrences of TC in the sequence."""
    logger.debug("Counting TC sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('TC', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_AC_sequences(sequence):
    """Count occurrences of AC in the sequence."""
    logger.debug("Counting AC sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('AC', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_CC_sequences(sequence):
    """Count occurrences of CC in the sequence."""
    logger.debug("Counting CC sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find('CC', start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_NC_sequences(sequence, context):
    """Count occurrences of the specified context (AC,CC,TC) in the sequence."""
    logger.debug(f"Counting {context} sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find(context, start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1
    return count

def count_GN_sequences(sequence, context):
    """Count occurrences of the specified context (GA, GT, GG) in the sequence."""
    logger.debug(f"Counting {context} sequences.")
    count = 0
    start_index = 0
    while True:
        index = sequence.find(context, start_index)
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

def find_N_positions(window, start_position, context):
    """Find positions of specified contexts (GA, GT, GG,TC,AC,CC) in the window relative to the mtDNA sequence."""
    logger.debug(f"Finding {context} contexts in the window.")
    positions = []
    start_index = 0
    while True:
        index = window.find(context, start_index)
        if index == -1:
            break
        positions.append(start_position + index + 1)  # +1 for 1-indexed
        start_index = index + 1
    return positions

def process_mtDNA(mtDNA_seq, pos):
    """Main function which processes the DNA"""
    logger.info("Processing mtDNA sequence for position %d.", pos)

    nospace_mtDNA = capitalize(remove_whitespace(mtDNA_seq))
    consecutive_TC_positions = find_consecutive_TC_sequences(nospace_mtDNA)
    consecutive_CC_positions = find_consecutive_CC_sequences(nospace_mtDNA)
    consecutive_AC_positions = find_consecutive_AC_sequences(nospace_mtDNA)
    consecutive_GA_positions = find_consecutive_GA_sequences(nospace_mtDNA)
    consecutive_GT_positions = find_consecutive_GT_sequences(nospace_mtDNA)
    consecutive_GG_positions = find_consecutive_GG_sequences(nospace_mtDNA)
    dummy = 0
    
    if pos in consecutive_TC_positions:
        logger.info("Base at position %d is in a 5'-TC context.", pos)
        ref, mut, all_windows, dum = 'C', 'T', [], []
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
        FLAG=None
        FLAG=None
        if right_adjacent_bases[0] == 'C':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='C':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            TC_windows = generate_NC_windows(circular_seq, pos, window_size)
            TALES = False
            for num, window in enumerate(TC_windows, start=4):
                window_description = f"Position {num} from the 3' end"
                ws = f"{window_size}bp"
                off_target_sites_topstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_bottomstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                #marked_window = mark_bases(window, window_size - num + 1, [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])])
                marked_window = mark_bases(
                        window, 
                        window_size - num + 1, 
                        [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(window[3:8])]
                )
                int_pos = marked_window.find(']')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos + 1] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos + 1)
                if int_pos - 3 >= 0 and marked_window[int_pos - 3] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos - 3)

                start_position = pos - (window_size - num)  # Adjust this according to your indexing logic
                tc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GG')
                # Combine positions into the respective lists
                if (pos - 1) in tc_positions or (pos - 1) in ac_positions or (pos - 1) in cc_positions or (pos + 1) in tc_positions or (pos + 1) in ac_positions or (pos + 1) in cc_positions:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                else:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy
                if pos in ftc:
                    ftc.remove(pos) #to remove the target position itself
                    
                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))
    
    elif pos in consecutive_AC_positions:
        logger.info("Base at position %d is in a 5'-AC context.", pos)
        ref, mut, all_windows, dum = 'C', 'T', [], []
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
        FLAG=None
        if right_adjacent_bases[0] == 'C':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='C':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            AC_windows = generate_NC_windows(circular_seq, pos, window_size)
            TALES = False
            for num, window in enumerate(AC_windows, start=4):
                window_description = f"Position {num} from the 3' end"
                ws = f"{window_size}bp"
                off_target_sites_topstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_bottomstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                #marked_window = mark_bases(window, window_size - num + 1, [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])])
                marked_window = mark_bases(
                        window, 
                        window_size - num + 1, 
                        [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(window[3:8])]
                )
                int_pos = marked_window.find(']')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos + 1] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos + 1)
                if int_pos - 3 >= 0 and marked_window[int_pos - 3] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos - 3)
                start_position = pos - (window_size - num)  # Adjust this according to your indexing logic
                tc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GG')
                # Combine positions into the respective lists
                if (pos - 1) in tc_positions or (pos - 1) in ac_positions or (pos - 1) in cc_positions or (pos + 1) in tc_positions or (pos + 1) in ac_positions or (pos + 1) in cc_positions:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                else:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy
                if pos in ftc:
                    ftc.remove(pos) #to remove the target position itself

                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))

    elif pos in consecutive_CC_positions:
        logger.info("Base at position %d is in a 5'-CC context.", pos)
        ref, mut, all_windows, dum = 'C', 'T', [], []
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
        FLAG=None
        if right_adjacent_bases[0] == 'C':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='C':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            CC_windows = generate_NC_windows(circular_seq, pos, window_size)
            TALES = False
            for num, window in enumerate(CC_windows, start=4):
                window_description = f"Position {num} from the 3' end"
                ws = f"{window_size}bp"
                off_target_sites_topstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_bottomstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                #marked_window = mark_bases(window, window_size - num + 1, [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])])
                marked_window = mark_bases(
                        window, 
                        window_size - num + 1, 
                        [(x + 3) for x in find_consecutive_GA_sequences(window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(window[3:8])]
                )
                int_pos = marked_window.find(']')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos + 1] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos + 1)
                if int_pos - 3 >= 0 and marked_window[int_pos - 3] == 'C':
                    final_window = mark_base_at_position(final_window, int_pos - 3)
                
                

                start_position = pos - (window_size - num)  # Adjust this according to your indexing logic
                tc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(window[3:8], start_position + 3 -1, 'GG')
                # Combine positions into the respective lists
                if (pos - 1) in tc_positions or (pos - 1) in ac_positions or (pos - 1) in cc_positions or (pos + 1) in tc_positions or (pos + 1) in ac_positions or (pos + 1) in cc_positions:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                    #off_target_sites = off
                else:
                    ftc = tc_positions + ac_positions + cc_positions + dum
                    fga = ga_positions + gt_positions + gg_positions 
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy
                    #off_target_sites = off
                if pos in ftc:
                    ftc.remove(pos) #to remove the target position itself

                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))

    elif pos in consecutive_GA_positions:
        logger.info("Base at position %d is in a 5'-GA context.", pos)
        ref, mut, all_windows, dum = 'G', 'A', [], []
        comple_nospace_mtDNA = complementing(nospace_mtDNA)
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", reverse_complement(left_adjacent_bases[::-1]), reverse_complement(right_adjacent_bases[::-1]))
        FLAG=None
        if right_adjacent_bases[0] == 'G':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='G':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            TALES = False
            GA_windows = generate_GN_windows(comple_nospace_mtDNA + comple_nospace_mtDNA, pos, window_size)
            for num, window in enumerate(GA_windows, start=4):
                window_description = f"Position {num} from the 5' end"
                ws = f"{window_size}bp"
                off_target_sites_bottomstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_topstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                reverse_window = complementing(window[::-1]) #to get top strand
                #marked_window = mark_bases(reverse_window, num, [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])])
                marked_window = mark_bases(
                        reverse_window, 
                        num, 
                        [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(reverse_window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(reverse_window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(reverse_window[3:8])]
                )
                int_pos = marked_window.find('[')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos - 1] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos - 1)
                if int_pos + 3 >= 0 and marked_window[int_pos + 3] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos + 3)

                start_position = pos - num + 1 
                tc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GG')
                #logger.info('%s , %s , %s, %s ', pos-1, ga_positions,gt_positions,gg_positions)
                if (pos - 1) in ga_positions or (pos - 1) in gt_positions or (pos - 1) in gg_positions or (pos+1) in ga_positions or (pos+1) in gt_positions or (pos+1) in gg_positions:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                else:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy
                if pos in fga:
                    fga.remove(pos)
            
                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))

    elif pos in consecutive_GG_positions:
        logger.info("Base at position %d is in a 5'-GG context.", pos)
        ref, mut, all_windows, dum = 'G', 'A', [], []
        comple_nospace_mtDNA = complementing(nospace_mtDNA)
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", reverse_complement(left_adjacent_bases[::-1]), reverse_complement(right_adjacent_bases[::-1]))
        FLAG=None
        if right_adjacent_bases[0] == 'G':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='G':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            TALES = False
            GG_windows = generate_GN_windows(comple_nospace_mtDNA + comple_nospace_mtDNA, pos, window_size)
            for num, window in enumerate(GG_windows, start=4):
                window_description = f"Position {num} from the 5' end"
                ws = f"{window_size}bp"
                off_target_sites_bottomstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_topstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                reverse_window = complementing(window[::-1]) #to get top strand
                #marked_window = mark_bases(reverse_window, num, [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])])
                marked_window = mark_bases(
                        reverse_window, 
                        num, 
                        [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(reverse_window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(reverse_window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(reverse_window[3:8])]
                )
                int_pos = marked_window.find('[')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos - 1] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos - 1)
                if int_pos + 3 >= 0 and marked_window[int_pos + 3] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos + 3)
                start_position = pos - num + 1 
                tc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GG')

                if (pos - 1) in ga_positions or (pos - 1) in gt_positions or (pos - 1) in gg_positions or (pos+1) in ga_positions or (pos+1) in gt_positions or (pos+1) in gg_positions:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                else:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy

                if pos in fga:
                    fga.remove(pos)
                
                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))

    elif pos in consecutive_GT_positions:
        logger.info("Base at position %d is in a 5'-GT context.", pos)
        ref, mut, all_windows, dum = 'G', 'A', [], []
        comple_nospace_mtDNA = complementing(nospace_mtDNA)
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - (16 + 15)
        end_index = pos + (15 + 15)
        adjacent_bases = circular_seq[start_index:end_index]
        marked_adjacent = mark_bases(adjacent_bases, 31, find_consecutive_GA_sequences(adjacent_bases) + find_consecutive_GG_sequences(adjacent_bases) + find_consecutive_GT_sequences(adjacent_bases) + find_consecutive_TC_sequences(adjacent_bases) + find_consecutive_AC_sequences(adjacent_bases) + find_consecutive_CC_sequences(adjacent_bases))
        #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", reverse_complement(left_adjacent_bases[::-1]), reverse_complement(right_adjacent_bases[::-1]))
        FLAG=None
        if right_adjacent_bases[0] == 'G':
            dummy = 1
            dum.append(pos+1)
            FLAG=True
        if left_adjacent_bases[-1] =='G':
            dummy += 1
            dum.append(pos-1)
            FLAG=True
        for window_size in range(14, 19):
            TALES = False
            GT_windows = generate_GN_windows(comple_nospace_mtDNA + comple_nospace_mtDNA, pos, window_size)
            for num, window in enumerate(GT_windows, start=4):
                window_description = f"Position {num} from the 5' end"
                ws = f"{window_size}bp"
                off_target_sites_bottomstrand = count_TC_sequences(window[-8:-3]) + count_CC_sequences(window[-8:-3])  + count_AC_sequences(window[-8:-3]) #finds the off-target sites in each window
                off_target_sites_topstrand = count_GA_sequences(window[3:8]) + count_GG_sequences(window[3:8]) + count_GT_sequences(window[3:8]) #finds the off-target sites in each window
                #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy #adding both the off_target sites from both the strands
                reverse_window = complementing(window[::-1]) #to get top strand
                #marked_window = mark_bases(reverse_window, num, [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])])
                marked_window = mark_bases(
                        reverse_window, 
                        num, 
                        [(x + 3) for x in find_consecutive_GA_sequences(reverse_window[3:8])] + 
                        [(x + window_size - 8) for x in find_consecutive_TC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_AC_sequences(reverse_window[-8:-3])] + 
                        [(x + window_size - 8) for x in find_consecutive_CC_sequences(reverse_window[-8:-3])] + 
                        [(x + 3) for x in find_consecutive_GG_sequences(reverse_window[3:8])] + 
                        [(x + 3) for x in find_consecutive_GT_sequences(reverse_window[3:8])]
                )
                int_pos = marked_window.find('[')
                final_window = marked_window  # Start with the original marked_window
                if int_pos + 1 < len(marked_window) and marked_window[int_pos - 1] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos - 1)
                if int_pos + 3 >= 0 and marked_window[int_pos + 3] == 'G':
                    final_window = mark_base_at_position(final_window, int_pos + 3)

                start_position = pos - num + 1 
                tc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'TC')
                ac_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'AC')
                cc_positions = find_N_positions(reverse_window[-8:-3], start_position + window_size - 8, 'CC')
                ga_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GA')
                gt_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GT')
                gg_positions = find_N_positions(reverse_window[3:8], start_position + 3 -1, 'GG')

                if (pos - 1) in ga_positions or (pos - 1) in gt_positions or (pos - 1) in gg_positions or (pos+1) in ga_positions or (pos+1) in gt_positions or (pos+1) in gg_positions:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand 
                else:
                    ftc = tc_positions + ac_positions + cc_positions 
                    fga = ga_positions + gt_positions + gg_positions + dum
                    #off_target_sites = off_target_sites_topstrand + off_target_sites_bottomstrand + dummy

                if pos in fga:
                    fga.remove(pos)
                off = final_window.count('{')
                sorted_positions = sorted(set(ftc + fga))  # Create a sorted version of the combined list
                all_windows.append((PIPELINE, pos, ref, mut, ws, final_window, window_description, off, sorted_positions, TALES, FLAG))
    else:
        logger.warning("Base at position %d is not in a editable context and cannot be edited by the %s pipeline.", pos, PIPELINE)
        print(f"Position {pos} is not in a editable context and cannot be edited by the {PIPELINE}.")
        return [], []  # Return empty lists to indicate failure

    return all_windows, adjacent_bases

def append_to_excel(all_windows, additional_file, output_file):
    """Search positions from 'ftc+fga' in annotated Excel file and append to all_windows."""
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
    parser.add_argument('additional_file', type=str, help='Excel file containing additional bystander information') #this can be hardcoded ig
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        logger.error("The file %s does not exist.", args.input_file)
        return
    
    if not os.path.isfile(args.additional_file):
        logger.error("The additional bystander file - %s does not exist.", args.additional_file)
        return
    
    logger.info("Reading the input sequence %s.", args.input_file)
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
