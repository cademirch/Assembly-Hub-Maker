# Assembly Hub Maker
A python script for making UCSC Genome Browser assembly hubs from NCBI data.  

This script takes the accession of a genome as input and downloads the requsite files for the assembly hub from the NCBI RefSeq database. The script also requires the path to a directory containing the file `genomes.txt`, which the script will edit with the assembly hub information.

# Usage:
`python assembly_hub.py -g [accession] -p [path to genomes.txt]`  

