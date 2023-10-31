# SeqAnalyzer

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://seq-analyzer.streamlit.app/)
[![python](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)

Website: https://seq-analyzer.streamlit.app/

SeqAnalyzer is a web application that allows you to analyze and compare sequences in three steps:

1. Provide sequences in [FASTA format](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/)
   - Example files from [NCBI](https://www.ncbi.nlm.nih.gov/) for Human coronavirus NL63 [GCF_000853865](example_files/GCF_000853865.fna) and SARS-CoV-2 [NC_045512](example_files/NC_045512.fasta) have been provided in the `example_files` directory of this repository
2. Select a function to be implemented on your sequence(s)
3. Click the "Analyze" button to obtain your results

# Functions

SeqAnalyzer currently supports the following functions.

## Transcription

The transcription function returns the mRNA sequence of the provided DNA sequence. The transcribed mRNA sequence can be downloaded as a text file.

## Translation

The translation function returns the protein sequence of the provided DNA sequence. By default, translation will use the [standard genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) from NCBI. The translated protein sequence can be downloaded as a text file.
