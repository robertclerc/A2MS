#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:30:16 2018

@author: clerc
"""

import os
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



#Objective

#Choose the refseq genes and transcrits for each consequence in mismatchs genome
#annotate the mismatchs with domains and genes names
#reverse the reads when the gene is complementary



def converse(nt):
    ''' documentation in preparation '''
    retour = []

    for n in nt.split(","):
        
        retour.append(str(Seq(n, IUPAC.unambiguous_dna).reverse_complement()))

    return ",".join(retour)

def choix_ref_seq(CSQ):
    ''' documentation in preparation '''
    dic_CSQ = {}

    for i in CSQ:

        info = i.split("|")
        gene = info[1]
        transcrit = info[2]

        if (transcrit in dic_CSQ) == False:

            dic_CSQ[transcrit] = {}
            dic_CSQ[transcrit]["refseq"] = 0
            dic_CSQ[transcrit]["info"] = ""


        dic_CSQ[transcrit]["refseq"] += dic_ref_seq[gene][transcrit]
        dic_CSQ[transcrit]["info"] = i


    max_ = 0
    max_transcrit = ""

    for transcrit in dic_CSQ:

        if dic_CSQ[transcrit]["refseq"] >= max_:

            max_ = dic_CSQ[transcrit]["refseq"]
            max_transcrit = dic_CSQ[transcrit]["info"]

    CSQ = [max_transcrit]

    return CSQ

def choix_CSQ_consequence(CSQ):
    ''' documentation in preparation '''
    liste_cons = []

    liste_CSQ = []

    retour = []


    for info in CSQ:

        cons = info.split("|")[7]

        if cons not in liste_cons:

            liste_cons.append(cons)

    if len(liste_cons) > 1:    

        for cons in liste_cons:

            liste_CSQ = []

            for info in CSQ:

                cons_info = info.split("|")[7]

                if cons == cons_info:

                    liste_CSQ.append(info)

            retour.extend(choix_ref_seq(liste_CSQ))

        CSQ = list(retour)

    else:

        CSQ = choix_ref_seq(CSQ)

    return CSQ

def file_ensembl(path):
    #ok for domaine
    #ok for transmembranaire
    #ok for gene_name
    ''' documentation in preparation '''
    dic_retour = {}
    dic_entete = {}

    list_file = os.listdir(path)

    for i in range(len(list_file)):

        files = open(path+list_file[i])

        read = files.readline()

        while read:

            row = read[:-1].split(",")

            if "Gene stable ID" in row[0] or "Ensembl Gene ID" in row[0]:

                #Building header

                for j in range(2,len(row)):

                    dic_entete[j] = row[j]
   
            else:

                gene = row[0]
                transcrit = row[1]

                if (gene in dic_retour) == False:

                    dic_retour[gene]={}

                if (transcrit in dic_retour[gene]) == False:

                    dic_retour[gene][transcrit] = []

                for j in dic_entete:

                    if row[j] != "":

                        dic_retour[gene][transcrit].append(row[j])

            read=files.readline()

        files.close()
        dic_entete={}
        
    return dic_retour

def file_refseq(path):
    ''' documentation in preparation '''
    #ok for refseq

    dic_retour = {}
    dic_entete = {}

    list_file = os.listdir(path)

    for i in range(len(list_file)):

        files = open(path+list_file[i])

        read = files.readline()

        while read:

            row = read[:-1].split(",")

            if "Gene stable ID" in row[0]:

                #Building header

                for j in range(2,len(row)):

                    dic_entete[j] = row[j]

            else:

                gene = row[0]
                transcrit = row[1]

                if (gene in dic_retour) == False:

                    dic_retour[gene] = {}

                if (transcrit in dic_retour[gene]) == False:

                    dic_retour[gene][transcrit] = []

                for j in dic_entete:

                    if row[j] !="":

                        dic_retour[gene][transcrit].append(row[j])

            read=files.readline()

        files.close()
        dic_entete={}

    for gene in dic_retour:

        for transcrit in dic_retour[gene]:

            dic_retour[gene][transcrit]=len(dic_retour[gene][transcrit])

    return dic_retour



parser = argparse.ArgumentParser()

choice = parser.add_argument_group('Choix cohorte')

choice.add_argument('-i','-path_in',type=str,help='Fichier d entree',required=True)
choice.add_argument('-o','-path_out',type=str,help='Fichier de sortie',required=True)
choice.add_argument('-e','-path_Ensembl',type=str,help='Chemin du repertoire Ensembl',required=True)

args = parser.parse_args()

path_e=args.e

path_transm=path_e+"transmembranaire/"

dic_transm=file_ensembl(path_transm)

path_ref_seq=path_e+"refseq/"

dic_ref_seq=file_refseq(path_ref_seq)

path_gene_name=path_e+"gene_name/"

dic_gene_name=file_ensembl(path_gene_name)

path_domain=path_e+"domaine/"

dic_domain=file_ensembl(path_domain)   

path_PDB=path_e+"PDB/"

dic_PDB=file_ensembl(path_PDB)   

path_in=args.i
path_out=args.o

path_in_enrichissement="/".join(args.o.split("/")[0:8])+"/"+"gene_name.txt"

f_in=open(path_in,"r")
f_out=open(path_out,"w")
f_out_2=open(path_in_enrichissement,"w")

f_out.write("##CHROMOSOME\tN_SNP\tN_gene\tInfo\tSNP_Donneur\tSNP_Receveur\tread_d\tread_r\ttype_nt_d\ttype_nt_r\tQual\tref_vcf\talt_vcf\tImput_d\tImput_r\n")
       
read=f_in.readline()

while read:

    if "#" not in read[0]:

        row=read[:-1].split("\t")
        chromosome=row[0]
        SNP=row[1]
        liste_CSQ=row[3].split(",")
        nt_donneur=row[4]
        nt_receveur=row[5]

        read_donneur=row[6]
        read_receveur=row[7]
        qual=row[8]

        type_nt_donneur=row[9]
        type_nt_receveur=row[10]

        new_nt_donneur=nt_donneur
        new_nt_receveur=nt_receveur
        
        ref_vcf=row[11]
        alt_vcf=row[12]
        
        Input_d=row[13]
        Input_r=row[14]

        

        #choose refseq

        if len(liste_CSQ)>1:

            liste_CSQ=choix_CSQ_consequence(liste_CSQ)

        for i in range(len(liste_CSQ)):
            
           

            gene=liste_CSQ[i].split("|")[1]
            transcrit=liste_CSQ[i].split("|")[2]
            sens=int(liste_CSQ[i].split("|")[9])
            CSQ=liste_CSQ[i]

            if (gene in dic_transm)==True:

                if (transcrit in dic_transm[gene])==True:

                    transm="transmembranaire="+".".join(list(set(dic_transm[gene][transcrit])))

                    if transm.split("=")[1]==[] or transm.split("=")[1]=="":

                        transm="transmembranaire=None_transm"
                    
                else:

                    transm="transmembranaire=None_transm"
                    
            else:

                transm="transmembranaire=None_transm"

            CSQ=CSQ+"|"+transm

            if (gene in dic_domain)==True:

                if (transcrit in dic_domain[gene])==True:

                    dom="domaine="+".".join(list(set(dic_domain[gene][transcrit])))

                    if dom.split("=")[1]==[]:

                        dom="domaine=None_domaine"

                else:

                    dom="domaine=None_domaine"

            else:

                dom="domaine=None_domaine"

            CSQ=CSQ+"|"+dom
 
            if (gene in dic_PDB)==True:

                if (transcrit in dic_PDB[gene])==True:

                    PDB="PDB="+".".join(list(set(dic_PDB[gene][transcrit])))

                    if PDB.split("=")[1]=="":

                        PDB="PDB=None_PDB"

                else:

                    PDB="PDB=None_PDB"

            else:

                PDB="PDB=None_PDB"

            CSQ=CSQ+"|"+PDB

            if (gene in dic_gene_name)==True:

                if (transcrit in dic_gene_name[gene])==True:
                    
                    gene_name="gene_name="+dic_gene_name[gene][transcrit][0]

                    if gene_name.split("=")[1]==[]:

                        gene_name="gene_name=None_name"

                else:

                    gene_name="gene_name=None_name"
            else:

                gene_name="gene_name=None_name"

            CSQ=CSQ+"|"+gene_name
            
            

            if sens==-1:
                


                if nt_donneur!="None":

                    new_nt_donneur=converse(nt_donneur)
                    

                if nt_receveur!="None":    

                    new_nt_receveur=converse(nt_receveur)
                    
                    
                ref_vcf=converse(ref_vcf)
                alt_vcf=converse(alt_vcf)

 
                
                    

            else:

                 new_nt_donneur=nt_donneur
                 new_nt_receveur=nt_receveur
                 

            f_out_2.write(gene_name+"\n")
            f_out.write(chromosome+"\t"+SNP+"\t"+gene_name+"\t"+CSQ+"\t"+new_nt_donneur+"\t"+new_nt_receveur+"\t"+read_donneur+"\t"+read_receveur+"\t"+type_nt_donneur+"\t"+type_nt_receveur+"\t"+qual+"\t"+ref_vcf+"\t"+alt_vcf+"\t"+str(Input_d)+"\t"+str(Input_r)+"\n")                                                           

    read=f_in.readline()

f_in.close()
f_out.close()
f_out_2.close()