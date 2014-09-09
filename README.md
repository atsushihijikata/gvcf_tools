gvcf_tools
==========

handling with gVCF (genome variant call format) files

Converting gVCF format into several formats.

Currently two formats are supported.
 1. multi-FASTA format with variations in the gVCF file.
 2. conventional VCF format.

Prerequisite:
    1. gVCF format file
    2. genome data (FASTA format) used for mapping
    
See gvcf_tools.py for details.

Usage:  
  python gvcf_tools.py sample.g.vcf

For more options:
  python gvcf_tools.py -help  
