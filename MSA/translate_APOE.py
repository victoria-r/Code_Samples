#!/usr/bin/env python3
# translate_APOE.py

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Get fasta file
fasta_infile = "APOE_refseq_transcript.fasta"

trim_aa = []

with open("apoe_aa.fasta", "w") as apoe_aa:
    # Parse fasta file
    for seq_record in SeqIO.parse(fasta_infile, "fasta"):
        
        # Translate sequences
        dna = seq_record.seq.translate()
            
        # Make SeqRecord object for writing
        protein = SeqRecord(dna, id=seq_record.id, description="Translated orthologs")
   
        # Write tranlsated sequence to fatsa
        SeqIO.write(protein, apoe_aa, "fasta")
