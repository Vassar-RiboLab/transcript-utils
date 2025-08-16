from Bio import SeqIO
import pandas as pd
import numpy as np
import regex
import re

def find_best_match(long_string, short_string, max_allowable_distance):
    """
    Finds the substring in the long_string with the smallest edit distance,
    up to a maximum allowable distance.

    Args:
        long_string (str): The string to search within.
        short_string (str): The string to find the best match for.
        max_allowable_distance (int): The maximum edit distance allowed for a match.

    Returns:
        A tuple containing the best matching substring, its start and end indices,
        and the edit distance. Returns None if no match is found within the
        specified maximum distance.
    """
    # Iterate from distance 0 up to the maximum allowed.
    for distance in range(max_allowable_distance + 1):
        # The pattern looks for a match with AT MOST 'distance' errors.
        pattern = f"({short_string}){{e<={distance}}}"
        match = regex.search(pattern, long_string)

        if match:
            # Because we start with distance=0, the first match found
            # is guaranteed to have the smallest possible edit distance.
            actual_distance = sum(match.fuzzy_counts)
            return (match.group(0), match.span(), actual_distance)

    # If the loop completes, no match was found within the allowed distance.
    return None

def has_overlap(main_string: str, sub1: str, sub2: str) -> bool:
    """
    Checks if two substrings overlap within a main string.

    This function assumes both substrings are present in the main string.
    """
    # Find the starting index of both substrings
    start1 = main_string.find(sub1)
    start2 = main_string.find(sub2)

    # Determine which substring appears first and assign start/end positions
    if start1 < start2:
        # sub1 comes first
        first_sub_end = start1 + len(sub1)
        second_sub_start = start2
    else:
        # sub2 comes first (or they start at the same index)
        first_sub_end = start2 + len(sub2)
        second_sub_start = start1

    # The core check: Does the first substring's end index surpass
    # the second substring's start index? If so, they overlap.
    return first_sub_end > second_sub_start


# Global var - length of window around start and stop codon to search for match
segment_len = 15

def extract_cds_info(transcript: str, cds_ids: list[str], cds_seqs: list[str], exact_list: list[bool]):
    """
    Removes all CDS regions within a transcript and returns a list
    of non-CDS elements in order, along with adjacency and overlap information.

    Handles contiguous and overlapping CDS's.

    :param transcript: String of nucleotides containing CDS's.
    :param cds_seqs: List of CDS sequences to be found in the transcript.
    :param exact_list: List of boolean values where True indicates exact match in transcript.
    :return: A tuple containing:
        - ordered_cds_ids (list[str]): CDS ids in order of appearance.
        - ordered_cds_seqs (list[str]): CDS's in order of appearance.
        - ordered_exact_list (list[bool]): True if exact match in transcript, order of appearance.
        - cds_len_diff (list[int]): Found CDS length minus catalogued CDS length.
        - cds_start_index (list[int]): Transcript indices of start codon ATG.
        - cds_stop_index (list[int]): Transcript indices of stop codon TAA, TAG, or TGA.
        - cds_in_frame (list[bool]): True if found cds length is multiple of 3.
        - non_cds_seqs (list[str]): All non-CDS regions in order. Empty strings for missing regions,
                                    in case of overlap, missing 5'UTR or 3'UTR or other similar phenomena.
        - contiguous_list (list[bool]): True if a CDS is contiguous to the next.
        - overlap_list (list[bool]): True if a CDS overlaps with the next.
        - frame_next_list (list[int]): Values from 0 to 2 to indicate frame of CDS with next.
        - fputr (str): The 5' Untranslated Region.
        - tputr (str): The 3' Untranslated Region.
    """
    transcript = transcript.upper()

    # Parameter for approximate search
    max_distance = 2

    # Return early if there are no CDS sequences to find.
    if not cds_seqs: 
        # The entire transcript is a non-CDS region.
        non_cds = [transcript] if transcript else []
        return ([], [], [], [], [], [], [], non_cds, [], [], [], "", "")

    # Step 1: Find the start and end of each CDS. Store as objects to keep data together.
    exact_indices = [i for i, x in enumerate(exact_list) if x]
    cds_info_list = []
    for i in range(len(cds_ids)):
        cds_id = cds_ids[i]
        cds = cds_seqs[i]
        cds_upper = cds.upper()

        start_seq = cds_upper[:segment_len]
        start_index = transcript.find(start_seq)

        end_seq = cds_upper[(len(cds) - segment_len):] 
        end_index = transcript.find(end_seq) + segment_len 

        if start_index == -1 and end_index == -1:
            match_info = find_best_match(transcript, cds, max_distance)
            if not match_info:
                continue
            _, (start_index,end_index), _ = match_info

        # In the event that there is a small difference in one of the segments, use find_best_match 
        # Sometimes there is long sequence interjected between start/stop codon and rest of cds
        # To account for this, get range where start/stop probably exists and find nearest start/stop
        extend=100

        if start_index == -1:
            start_match = find_best_match(transcript, start_seq, 1)
            if not start_match:
                start_codon = "ATG"

                # Extend before theoretical start index (based on CDS len) by certain num. of nucleotides
                theor_start_index = end_index - len(cds)
                search_seq = transcript[theor_start_index-extend:theor_start_index]

                start_index = search_seq.rfind(start_codon)

            else: _, (start_index, _), _ = start_match

        if end_index == -1:
            end_match = find_best_match(transcript, end_seq, 1)
            if not end_match:
                stop_codon = cds[len(cds)-3:]

                # Extend after theoretical end index (based on CDS len) by certain num. of nucleotides
                theor_end_index = start_index + len(cds)
                search_seq = transcript[theor_end_index:theor_end_index+extend]

                end_index = search_seq.find(stop_codon)+3

                if end_index==-1:
                    continue
            else: _, (_, end_index), _ = start_match

        # Extract the difference in lengths between the found cds and the catalogued cds
        match_len = end_index - start_index 
        diff = match_len - len(cds)

        # Check that the found cds is in frame
        # Future implementation: If not, want to find next longest ORF
        in_frame = match_len % 3 == 0

        exact = exact_list[i]

        if i not in exact_indices:
            for j in exact_indices:
                exact_cds = cds_seqs[j]
                exact_cds_upper = exact_cds.upper()
                exact_start = transcript.find(exact_cds_upper[:segment_len])
                exact_end = transcript.find(exact_cds_upper[(len(exact_cds) - segment_len):]) + segment_len

                if exact_start == start_index and exact_end == end_index:
                    continue

        
        found_cds = transcript[start_index:end_index+1]

        if not abs(diff) > len(cds*2):
            cds_info_list.append({
                "id": cds_id,
                "recorded_seq": cds,
                "found_seq": found_cds,
                "start": start_index,
                "end": end_index,
                "diff": diff,
                "in_frame": in_frame,
                "exact": exact
            })

    # If no provided CDS sequences were found in the transcript.
    if not cds_info_list:
        non_cds = [transcript] if transcript else []
        return ([], [], [], [], [], [], [], non_cds, [], [], [], "", "")

    # Step 2: Sort the CDSs based on their starting position.
    ordered_cds_info = sorted(cds_info_list, key=lambda x: x["start"])
    ordered_cds_ids = [info["id"] for info in ordered_cds_info]
    ordered_cds_seqs = [info["found_seq"] for info in ordered_cds_info]
    ordered_exact_list = [info["exact"] for info in ordered_cds_info]

    # Initialize return values
    non_cds_seqs = []
    contiguous_list = []
    overlap_list = []
    frame_next_list = []

    # Step 3: Extract the 5' UTR (the sequence before the first CDS).
    first_cds_start = ordered_cds_info[0]["start"]
    fputr = transcript[:first_cds_start]
    if fputr:
        non_cds_seqs.append(fputr)
    else: 
        non_cds_seqs.append("")

    # Step 4: Iterate through ordered CDS pairs to find introns and check for adjacency/overlap.
    # The pointer `last_cds_end` tracks our position after each processed CDS.
    last_cds_end = ordered_cds_info[0]["end"]

    for i in range(len(ordered_cds_info) - 1):
        current_cds_info = ordered_cds_info[i]
        next_cds_info = ordered_cds_info[i+1]

        # Check for contiguousness (no gap, no overlap)
        is_contiguous = current_cds_info["end"] == next_cds_info["start"]
        contiguous_list.append(is_contiguous)

        # Check for overlap
        is_overlap = current_cds_info["end"] > next_cds_info["start"]
        overlap_list.append(is_overlap)

        # Check the frame.
        frame = np.abs(next_cds_info["start"] - current_cds_info["end"]) % 3
        frame_next_list.append(frame)

        # An intron is the non-CDS region between two CDSs. It only exists if there is a gap.
        if current_cds_info["end"] < next_cds_info["start"]:
            intron = transcript[current_cds_info["end"]:next_cds_info["start"]]
            non_cds_seqs.append(intron)
        else:
            non_cds_seqs.append("")
        
        # To handle cases where CDSs are nested or have complex overlaps,
        # we advance our pointer to the furthest point reached so far.
        last_cds_end = max(last_cds_end, next_cds_info["end"])

    cds_len_diff = [info["diff"] for info in ordered_cds_info]
    cds_start_index = [info["start"] for info in ordered_cds_info]
    cds_stop_index = [info["end"] - 3 for info in ordered_cds_info]
    cds_in_frame = [info["in_frame"] for info in ordered_cds_info]

    # Step 5: Extract the 3' UTR (the sequence after the last point covered by a CDS).
    tputr = transcript[last_cds_end:]
    if tputr:
        non_cds_seqs.append(tputr)
    else:
        non_cds_seqs.append("")
    
    
    return (ordered_cds_ids, ordered_cds_seqs, ordered_exact_list, cds_len_diff, cds_start_index, cds_stop_index,
            cds_in_frame, non_cds_seqs, contiguous_list, overlap_list, frame_next_list, fputr, tputr)


def find_CDS(transcriptome: str, cds_fasta: str):
    """For each transcript in wf-transcriptomes generated transcriptome, finds each cds present and transcript elements not present in cds sequences.
    The purpose of this is to understand how the transcriptome is created, whether we can gather useful information from it
    in a de novo analysis, what is in the MSTRG sequences, and how alignment might be informed or misinformed by the generated transcriptome.

    Assumes that there is only one instance of each cds per transcript.

    :param transcriptome: Location of transcriptome fasta to be iterated through.
    :param cds_fasta: Location of NCBI CDS fasta - file containing all relevant CDS's and their ids.
                      For now only works when IDs have identifier "locus_tag=(CDS_ID)."
    :return: dataset of transcriptome match info
             Includes: {id(rowname):str, transcript:str, num_cds:int, num_cds_exact:int, 
             exact_list:[bool], cds_ids:[str], cds_seqs:[str], cds_len_diff:[int], cds_start_index:[int], 
             cds_stop_index:[int], cds_in_fram_list:[bool], non_cds_elems:[str], num_non_cds_elems:[str],
             contiguous_cds_list:[bool], overlap_cds_list:[bool], frame_next_cds_list:[int],
             fputr:str, tputr:str}
    """

    dict = {
        "id": [],
        "sequence": [],
        "num_cds": [],
        "num_cds_exact": [],
        "exact_list": [],
        "cds_ids": [],
        "cds_seqs": [],
        "cds_len_diff": [],
        "cds_start_index": [],
        "cds_stop_index": [],
        "cds_in_frame_list": [],
        "non_cds_elems": [],
        "num_non_cds_elems": [],
        "contiguous_cds_list": [],
        "overlap_cds_list": [],
        "frame_next_cds_list": [],
        "fputr": [],
        "tputr": [],
    }

    for T_record in SeqIO.parse(transcriptome, "fasta"):

        id = T_record.id

        transcript_seq = str(T_record.seq).upper()

        # Iterate over cds to find all cds included and non-cds elements in transcript
        # Can use the .replace method to get non-cds elements
        num_cds_exact = 0
        cds_ids = []
        cds_seqs = []
        exact_list = []
        for cds_record in SeqIO.parse(cds_fasta, "fasta"):

            cds_seq = str(cds_record.seq).upper()

            # Some cds's in transcripts have small insertions. Get first and last n nucleotides, 
            # if there is an exact match, defer to that. Otherwise, go with inexact
            if cds_seq in transcript_seq:
                # We have an exact match
                num_cds_exact += 1

                # No easy way to get the 'locus_tag', so we have to use re.search
                # cds_id = re.search(r"gene:(.*?)", cds_record.description).group(1)
                cds_id = re.search(r"\[locus_tag=(.*?)\]", cds_record.description).group(1)

                cds_ids.append(cds_id)
                cds_seqs.append(cds_seq)

                exact_list.append(True)

            elif cds_seq[:segment_len] in transcript_seq or cds_seq[len(cds_seq) - segment_len:] in transcript_seq:
                cds_id = re.search(r"\[locus_tag=(.*?)\]", cds_record.description).group(1)

                cds_ids.append(cds_id)
                cds_seqs.append(cds_seq)

                exact_list.append(False)
        
        # Use extract_cds_infor to collect non-CDS elements accounting for overlap and direct adjacency
        (ordered_cds_ids,
         ordered_cds_seqs,
         ordered_exact_list,
         cds_len_diff,
         cds_start_index,
         cds_stop_index,
         cds_in_frame_list,
         non_cds_seqs_inclusive, # Includes empty CDS regions from phenomena like overlaps
         contiguous_list,
         overlap_list,
         frame_list,
         fputr,
         tputr) = extract_cds_info(transcript_seq, cds_ids, cds_seqs, exact_list)

        num_cds = len(ordered_cds_ids)

        print(f"{num_cds} CDS's found in transcript {id}.")

        # Need to get non-inclusive non-CDS list to get actual number non-CDS
        non_cds_seqs = [seq for seq in non_cds_seqs_inclusive if seq != ""]

        # Now we have everything. Combine into the dictionary.
        dict["id"].append(id)
        dict["sequence"].append(transcript_seq)
        dict["num_cds"].append(len(ordered_cds_ids))
        dict["num_cds_exact"].append(num_cds_exact)
        dict["exact_list"].append(ordered_exact_list)
        dict["cds_ids"].append(ordered_cds_ids)
        dict["cds_seqs"].append(ordered_cds_seqs)
        dict["cds_len_diff"].append(cds_len_diff)
        dict["cds_start_index"].append(cds_start_index)
        dict["cds_stop_index"].append(cds_stop_index)
        dict["cds_in_frame_list"].append(cds_in_frame_list)
        dict["non_cds_elems"].append(non_cds_seqs_inclusive)
        dict["num_non_cds_elems"].append(len(non_cds_seqs))
        dict["contiguous_cds_list"].append(contiguous_list)
        dict["overlap_cds_list"].append(overlap_list)
        dict["frame_next_cds_list"].append(frame_list)
        dict["fputr"].append(fputr)
        dict["tputr"].append(tputr)
    
    # We have constructed the dictionary. Now, we make a dataframe.
    df = pd.DataFrame(dict)

    df = df.set_index('id')

    return df