gvcf_tools
==========

Handling with gVCF (genome variant call format) files
to convert gVCF into the other formats.

Currently two formats are available:

 1. multi-FASTA format with variations in the gVCF file.
 2. conventional VCF format.

Prerequisite:
    1. gVCF format file
    2. genome data (FASTA format) used for mapping
    
See gvcf_tools.py for details.

Usage:  
  type `python gvcf_tools.py sample.g.vcf`

For more options:  
  type `python gvcf_tools.py -help` 
