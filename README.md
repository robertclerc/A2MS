# A2MS Pipeline
<em> This software extract and build peptides around SNP differents between two patient from a WES 
represented by a Variant Calling Format (VCF) file, aligned with reference GRCh37. The goal is to product a fasta file with all differents peptides for further analysis.</em>

## Step 1 : compute_AMS_1.py

<em>According a list of SNP enable to create differences in protein sequences between to individuals called mismatch, this script detect all mismatch in the genome and extract relevant information for furthers analysis. This script reduced by about 90% the size of input files.</em>

## Step 2 : filter_and_annotations_2.py

<em>This script annote all mismatch detected ( refseq, transmembrane protein )</em>

## Step 3 : contruct_peptides_3.py

<em>Thanks to functions in fonctions_constructions.py, this script build the genetic sequences around mismatchs according the type of mismatch ( missense_variant, stop_deletion, etc ) and translate the genetic sequence to protein sequence for both individuals. See Donor_result.fasta and Recipient_result.fasta.</em>
