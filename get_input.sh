#!/bin/bash
# https://github.com/michalbukowski/tas-drug
# (c) 2023-2025 Michał Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# Department of Analytical Biochemistry, Faculty of Biochemistry, Biophysics and Biotechnology,
# Jagiellonian University, Krakow, Poland
# Distributed under GPL-3.0 license
# This notebook is a part of the TAS-Drug repository and is associated with the following study:
# Bukowski M, Banasik M, Chlebicka K, Bednarczyk K, Bonar E, Sokołowska D, Żądło T,
# Dubin G, Władyka B. Analysis of co-occurrence of type II toxin–antitoxin systems
# and antibiotic resistance determinants in Staphylococcus aureus. mSystems 0:e00957-24.
# https://doi.org/10.1128/msystems.00957-24

set -euo pipefail

mkdir -p input

echo "Downloading files to the 'input' directory..."

wget -q --show-progress -O input/genomes_list.txt 'https://ujchmura-my.sharepoint.com/:t:/g/personal/m_bukowski_uj_edu_pl/ER59axff3ydFhNb60uFRA2sBPrcYIyu6nYzQYxgmhMdgxw?e=rSVVjP&download=1'
wget -q --show-progress -O input/antitoxins.faa   'https://ujchmura-my.sharepoint.com/:u:/g/personal/m_bukowski_uj_edu_pl/ERJpKVgD7h9Mnb1hdMf0ymoB2ocydgRrqB9GKxuJJ1mUgg?e=icw6Na&download=1'
wget -q --show-progress -O input/toxins.faa       'https://ujchmura-my.sharepoint.com/:u:/g/personal/m_bukowski_uj_edu_pl/EU8INwBuPpFFtawYTfGai5UBuZ2tRKN6eaZsjlegpRKKXQ?e=a3PLqx&download=1'
wget -q --show-progress -O input/card.faa         'https://ujchmura-my.sharepoint.com/:u:/g/personal/m_bukowski_uj_edu_pl/EUHxUiL_Hb1JroaBF1qV9mIBkelcDkmCnpZd8FqxIY4VCw?e=mpqMjb&download=1'
wget -q --show-progress -O input/tblastn.tsv      'https://ujchmura-my.sharepoint.com/:u:/g/personal/m_bukowski_uj_edu_pl/EZoDlGsAITZMqJIKFcydeU4BK4AJAZflbiLBWeHHRSUcYw?e=kkFON3&download=1'

echo "Downloading successfully completed!"
