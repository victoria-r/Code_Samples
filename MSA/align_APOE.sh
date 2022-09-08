#!/usr/bin/env bash
# align_APOE.sh

# This script aligns the transcribed APOE orthologs using Clustal Omega

clustalo -i apoe_aa.fasta -t protein  -o aligned_apoe.fasta

