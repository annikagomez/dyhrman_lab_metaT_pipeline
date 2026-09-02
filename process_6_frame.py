#!/usr/bin/env python3

#Author: Annika Gomez, algomez@ldeo.columbia.edu
#May 18th, 2026

"""
Usage: The purpose of this script is to process the output of emboss transeq 
(full 6-frame translation of assembled transcripts) for downstream annotation.
Input: Full 6-frame translation in fasta format with one entry per frame (6* # transcripts)
Positional argument:
    1.  6-frame translation file name. If this is a path, the folder containing the 6-frame translation file 
    will be treated as the working directory
Optional argument:
    --min_len   Minimum peptide length in amino acids. Default = 50 
    --prefix    Prefix for output file names. Default = six_frame
Output: 
    1. Fasta file with one entry per transcript, all CDS longer than min_len cutoff concatenated and separated by *
    2. Fasta file with one entry per CDS, ID contains frame and coordinates in frame
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse


#function to separate coding sequences (split by '*') for each sequence
def get_CDS(seq,header,minimum):
    long_cds = []
    long_cds_seq = []
    CDS = seq.split('*')
    for i in CDS:
        if len(i) > minimum:
            long_cds.append(i)
            record = SeqRecord(Seq(i), id=header + '_' +str(seq.find(i)*3)+ '-' +str(seq.find(i)*3 + len(i)*3))
            long_cds_seq.append(record)
    return(long_cds, long_cds_seq)

#function to get combined and separate 6-frame translation for each transcript
#frames = list of each frame translation; trx = transcript name
def translate_trx(frames, trx,min_len):
    concat_trans = []
    all_records = []
    count = 1
    for frame in frames:
        seqs, seq_recs = get_CDS(frame, trx+'_'+str(count),min_len)
        count += 1
        all_records += seq_recs
        concat_trans += seqs
    all_trans = "*".join(concat_trans)
    if len(all_trans) > 1:
        one_record = SeqRecord(Seq(all_trans), id=trx)
    else:
        one_record = None
    return(all_records,one_record)

def process_6_frame(file,min_len,prefix):
    current_trx_frames = []
    all_cds = [] #individual coding regions longer than cutoff
    all_6_frame = [] #concatenated 6-frame translation of cds longer than cutoff
    if '/' in file:
        working_dir = file.rsplit('/',maxsplit=1)[0] + '/'
    else:
        working_dir = ''

    for index, record in enumerate(SeqIO.parse(file, "fasta")):
        current_trx_frames.append(str(record.seq))
        if (index+1)%6 == 0:
            all_frames, concat = translate_trx(current_trx_frames,record.id.rsplit('_',maxsplit=1)[0],min_len)
            all_cds += all_frames
            if isinstance(concat, SeqRecord):
                all_6_frame.append(concat)
            else:
                pass
            current_trx_frames = []
            if (index+1)%1000 == 0:
                print((index+1),flush=True)

    with open(working_dir+prefix+"_merged_cds.pep", "w") as output_handle:
        SeqIO.write(all_6_frame, output_handle, "fasta")

    with open(working_dir+prefix+"_separate_cds.pep", "w") as output_handle:
        SeqIO.write(all_cds, output_handle, "fasta")
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="6-frame processor", description="The purpose of this script is to process the output of emboss transeq \
        (full 6-frame translation of assembled transcripts) for downstream annotation.")
    parser.add_argument('filename')
    parser.add_argument('-l', '--min_len', dest='min_len', type=int, default = 50, help='Minimum peptide length in amino acids, default = 50')
    parser.add_argument('-p', '--prefix', dest='prefix', type=str, default = 'six_frame', help='Prefix for output files, default = six_frame')
    args = parser.parse_args()

    process_6_frame(args.filename,args.min_len,args.prefix)
