#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 13:25:23 2018

@author: clerc
"""


###########
#Objective#
###########

########################################################################
#prepare a pipeline to call all commands and script for compute AMS and#
############## extract the peptides around Mismatchs ###################
########################################################################

import argparse
import os


def open_init_file(path_init, list_all_pathv):
    """
    For prepare all paths variables
    """

    liste_retour = []

    dicc = locals()

    file_init = open(path_init, 'r')

    read = file_init.readline()

    while read:

        if "#" not in read[0]:

            name = read[:-1].split("=")[0]
            path = read[:-1].split("=")[1]

            dicc[name] = path

        read = file_init.readline()

    file_init.close()

    for path in list_all_pathv:

        liste_retour.append(dicc[path])
        

    return liste_retour

##########################
#tar.gz expected in input#
##########################

PARSER = argparse.ArgumentParser()

CHOICE = PARSER.add_argument_group('Choix cohorte')

CHOICE.add_argument('-d', '-path_in', type=str, help='File tar.gz of donor', required=True)

CHOICE.add_argument('-r', '-receveur', type=str, help='File tar.gz of recipient', required=True)

CHOICE.add_argument('-c', '-caller', type=str, choices=["HC", "platypus", "sambamba", "UG3", "merged"], help='Caller use', required=False, default='HC')

CHOICE.add_argument('-s_d', '-sexe_donneur', choices=["M", "F"], help='donor sex', required=True)

CHOICE.add_argument('-s_r', '-sexe_receveur', choices=["M", "F"], help='recipient sex', required=True)

CHOICE.add_argument('-i', '-init', type=str, help='Init files containing all paths of projects', required=True)

CHOICE.add_argument('-l', '-liste_lg', type=str, help='liste des longeures de peptides a construire', required=True)

CHOICE.add_argument('-e_d', '-ethnie_donneur', type=str, help='ethnie du donneur', required=True)

CHOICE.add_argument('-e_r', '-ethnie_receveur', type=str, help='ethnie du receveur', required=True)

CHOICE.add_argument('-CSQ', '-fields', type=str, choices=["CSQ", "VariantEffectPrediction"], help='sexe du receveur', required=False,default='CSQ')

CHOICE.add_argument('-rd', '-ratio_doneur', type=str , required=False,default=0.5)

CHOICE.add_argument('-rr', '-ratio_receveur', type=str, required=False,default=0.5)

ARGS = PARSER.parse_args()

#peptide size
LISTE_LONGUEURE = list(map(int, ARGS.l.split(",")))


#list of all variables paths
LIST_ALL_PATH = ["PATH_INPUT_D", "PATH_INPUT_R", "PATH_OUTPUT", "PATH_PROJET", "PATH_E", "PATH_REFSEQ","PATH_EXAC","PATH_GTEX","PATH_COHORTE","PATH_BED"]

PATH_INPUT_D, PATH_INPUT_R, PATH_OUTPUT, PATH_PROJET, PATH_E, PATH_REFSEQ, PATH_EXAC, PATH_GTEX, PATH_COHORTE, PATH_BED = open_init_file(ARGS.i, LIST_ALL_PATH)

PATH_PROJET="/home/clerc/Projet/A2MS/"


if ARGS.rd:
    
    ratio_donneur=ARGS.rd
    
    
if ARGS.rr:
    
    ratio_receveur=ARGS.rr

#Caller use
if ARGS.c:

    CALL = ARGS.c
    
ethnie_d=ARGS.e_d
ethnie_r=ARGS.e_r

SEXE_D = ARGS.s_d
SEXE_R = ARGS.s_r

#dir of patient
DIR_D=ARGS.d.split(".")[0]+"/"
DIR_R=ARGS.r.split(".")[0]+"/"


#patient's dir
FILE_D = PATH_INPUT_D+ARGS.d
FILE_R = PATH_INPUT_R+ARGS.r

#patient's name for path file
NAME_D = ARGS.d
NAME_R = ARGS.r



NAME_SIMPLE_D = NAME_D.split("_")[2].split(".")[0]
NAME_SIMPLE_R = NAME_R.split("_")[2].split(".")[0]


#merged name for path file
MERGED_NAME = NAME_SIMPLE_D+"_"+NAME_SIMPLE_R


PATH_OUTPUT_TABLE=PATH_COHORTE+"output_table/"+MERGED_NAME


if os.path.isdir(PATH_OUTPUT_TABLE)==False:
    
    os.system("mkdir "+PATH_OUTPUT_TABLE)

######################
#mkdir of dirs output#
######################

os.system("mkdir "+PATH_OUTPUT+MERGED_NAME)
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/temp")
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/plot_qual")
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/plot_qual/files")
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/plot_qual/plot")
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/table_premhc")
os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/out_MHC")

os.system("mkdir "+PATH_OUTPUT+MERGED_NAME+"/peptides")

os.system("mkdir "+PATH_COHORTE+"/output_table_premhc/"+MERGED_NAME)



#############
#gvcf to vcf#
#############

#os.system("bcftools convert --gvcf2vcf "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH.HC.clean_DP_and_chr.annot.g.vcf --fasta-ref /data2/clerc/Projet/Second_stack_cohort/grch37-reference-genome-ensembl-v1.fa --output-type z > "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_convertgvcf2vcf.vcf")
#os.system("bcftools convert --gvcf2vcf "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH.HC.clean_DP_and_chr.annot.g.vcf --fasta-ref /data2/clerc/Projet/Second_stack_cohort/grch37-reference-genome-ensembl-v1.fa --output-type z > "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_convertgvcf2vcf.vcf")


###############################
#extraction and gunzip of files#
################################

#if os.path.isfile(PATH_INPUT_D+ARGS.d):
#
#    os.system("gunzip "+PATH_INPUT_D+ARGS.d)
#
#else:
#
#    print("file no found or already decompressed")
#
#if os.path.isfile(PATH_INPUT_R+ARGS.r):
#
#    os.system("gunzip "+PATH_INPUT_R+ARGS.r)
#
#else:
#
#    print("file no found or already decompressed")


#################################
#Decompression des fichiers gvcf#
#################################

#os.system("gunzip "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH."+ARGS.c+".annot.g.vcf.gz")
#os.system("gunzip "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH."+ARGS.c+".annot.g.vcf.gz")



###############################
#Bed Tools intersect bed files#
###############################

#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_D)
#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_R)


#os.system("bedtools intersect -header -a /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_D+"/"+NAME_SIMPLE_D+"_hs37d5_BOTH.HC.annot.g.vcf -b /data2/clerc/Projet/Second_stack_cohort/bed_file/SureselectV5_Full_match.bed > /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_D+"/"+NAME_SIMPLE_D+"_hs37d5_BOTH.HC.annot_clean_DP_intersected_bed.g.vcf")
#os.system("bedtools intersect -header -a /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_R+"/"+NAME_SIMPLE_R+"_hs37d5_BOTH.HC.annot.g.vcf -b /data2/clerc/Projet/Second_stack_cohort/bed_file/SureselectV5_Full_match.bed > /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_R+"/"+NAME_SIMPLE_R+"_hs37d5_BOTH.HC.annot_clean_DP_intersected_bed.g.vcf")

#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_D)
#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_R)



#os.system("bedtools intersect -header -a /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_D+"/"+NAME_SIMPLE_D+"_convertgvcf2vcf.vcf-b /data2/clerc/Projet/Second_stack_cohort/bed_file/SureselectV5_Full_match.bed > /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_D+"/"+NAME_SIMPLE_D+"_convertgvcf2_intersected_bed.g.vcf")
#os.system("bedtools intersect -header -a /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_R+"/"+NAME_SIMPLE_R+"_convertgvcf2vcf.vcf -b /data2/clerc/Projet/Second_stack_cohort/bed_file/SureselectV5_Full_match.bed > /data2/clerc/Projet/Second_stack_cohort/vcf_input/deliverable_VARIANT_"+NAME_SIMPLE_R+"/"+NAME_SIMPLE_R+"_convertgvcf2_intersected_bed.g.vcf")

#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_D)
#os.system("python /data2/clerc/Projet/Second_stack_cohort/clean_gintersect_bed.py -p deliverable_VARIANT_"+NAME_SIMPLE_R)



##############
#tar of files#
##############
    

#os.system("tar xzf "+PATH_INPUT_D+ARGS.d[:-3]+" -C "+PATH_INPUT_D)
#os.system("tar xzf "+PATH_INPUT_R+ARGS.r[:-3]+" -C "+PATH_INPUT_R)
##

####bgzip of files
#os.system("bgzip "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH."+ARGS.c+".annot.vcf")
#os.system("bgzip "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH."+ARGS.c+".annot.vcf")
##
###tabix of files
#os.system("tabix -p vcf "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH."+ARGS.c+".annot.vcf.gz")
#os.system("tabix -p vcf "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH."+ARGS.c+".annot.vcf.gz")

#bgzip of files
#os.system("bgzip -f "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH."+ARGS.c+".annot_clean_DP_intersected_bed.g.vcf")
#os.system("bgzip -f "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH."+ARGS.c+".annot_clean_DP_intersected_bed.g.vcf")
###
####tabix of files
#os.system("tabix -f -p vcf "+PATH_INPUT_D+DIR_D+NAME_SIMPLE_D+"_hs37d5_BOTH.HC.annot_clean_DP_intersected_bed.g.vcf.gz")
#os.system("tabix -f -p vcf "+PATH_INPUT_R+DIR_R+NAME_SIMPLE_R+"_hs37d5_BOTH.HC.annot_clean_DP_intersected_bed.g.vcf.gz")



####################
#merge of vcf files#
####################
    

#os.system("vcf-merge " + PATH_INPUT_D+DIR_D+"/"+NAME_SIMPLE_D+"_hs37d5_BOTH."+ARGS.c+".annot_clean_DP_intersected_bed.g.vcf.gz " +PATH_INPUT_R+DIR_R+"/"+NAME_SIMPLE_R+"_hs37d5_BOTH."+ ARGS.c+".annot_clean_DP_intersected_bed.g.vcf.gz  | bgzip -c > "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.g.vcf.gz")


##decompression of merged file
#os.system("rm -r "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.vcf")

#os.system("gunzip "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.g.vcf.gz")

#os.system("gzip "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.g.vcf")

####################################
#annotation of merged file with vep#
####################################


#classic
#os.system("/data2/clerc/vep_dir/ensembl-vep/vep -i "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.g.vcf.gz --format vcf -o "+PATH_OUTPUT+MERGED_NAME+"/temp/out_annoted.g.vcf.gz --keep_csq --offline --force_overwrite --compress_output bgzip --assembly GRCh37 --vcf CSQ --fields Uploaded_variation,Gene,Feature,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,STRAND")

#print("/data2/clerc/vep_dir/ensembl-vep/vep -i "+PATH_OUTPUT+MERGED_NAME+"/temp/out_merged.g.vcf.gz --format vcf -o "+PATH_OUTPUT+MERGED_NAME+"/temp/out_annoted.g.vcf.gz --keep_csq --offline --force_overwrite --compress_output bgzip --assembly GRCh37 --vcf CSQ --fields Uploaded_variation,Gene,Feature,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,STRAND")




#####################
#Extract merged file#
#####################

#os.system("gunzip "+PATH_COHORTE+"vcf_output/"+MERGED_NAME+"/temp/out_annoted_MES_only.vcf.gz")


#############
#Compute AMS#
#############


#Compute AMS
os.system("python "+PATH_PROJET+"compute_AMS_1.py -c -v "+PATH_OUTPUT+MERGED_NAME+ "/temp/out_annoted.vcf.gz -d "+NAME_SIMPLE_D+" -r "+NAME_SIMPLE_R+ " -o "+PATH_OUTPUT +" -s_d " +SEXE_D +" -s_r "+SEXE_R +" -CSQ CSQ -pe "+PATH_E + " -min "+ str(10) +" -max "+ str(500) +" -PL "+ str(0.5) +" -rd " +str(ratio_donneur)+" -rr "+ str(ratio_receveur))


#filtering and annotation of SNPs

#####mistmatch file
os.system("python "+PATH_PROJET+"filter_and_annotation_2.py -i "+PATH_OUTPUT+MERGED_NAME+  "/temp/info_nt.txt -o "+PATH_OUTPUT+MERGED_NAME+"/temp/info_nt_filter.txt -e "+PATH_E)


#build peptides
os.system("python "+PATH_PROJET+"construct_peptides_3.py -e "+PATH_E+" -i "+ PATH_OUTPUT+MERGED_NAME+ " -l "+",".join(map(str, LISTE_LONGUEURE)) +" -e_d "+ethnie_d+" -e_r "+ethnie_r+" -e "+PATH_E + " -gtx "+ PATH_GTEX +" -p "+PATH_OUTPUT+MERGED_NAME +" -ex "+PATH_EXAC+" -t "+PATH_COHORTE+"output_table_premhc/"+MERGED_NAME + " -b "+PATH_BED)




