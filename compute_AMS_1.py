#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:47:24 2018

@author: clerc
"""



import argparse
import csv
import vcf
import os



#liste_all_consequence=['downstream_gene_variant', 'intron_variant', 'non_coding_transcript_variant', 'non_coding_transcript_exon_variant', 'splice_region_variant', 'upstream_gene_variant', 'synonymous_variant', 'missense_variant', '3_prime_UTR_variant', 'intergenic_variant', 'NMD_transcript_variant', '5_prime_UTR_variant', 'start_lost', 'inframe_insertion', 'inframe_deletion', 'splice_acceptor_variant', 'coding_sequence_variant', 'stop_gained', 'frameshift_variant', 'splice_donor_variant', 'mature_miRNA_variant', 'stop_lost', 'stop_retained_variant', 'start_retained_variant', 'protein_altering_variant', 'incomplete_terminal_codon_variant']

#liste_codant_SNP=['synonymous_variant', 'missense_variant', 'splice_region_variant&synonymous_variant', 'synonymous_variant&NMD_transcript_variant', 'missense_variant&NMD_transcript_variant', 'frameshift_variant', 'start_lost&NMD_transcript_variant', 'start_lost', 'missense_variant&splice_region_variant', 'inframe_deletion', 'splice_region_variant&synonymous_variant&NMD_transcript_variant', 'splice_acceptor_variant&coding_sequence_variant&intron_variant', 'stop_gained', 'inframe_insertion', 'missense_variant&splice_region_variant&NMD_transcript_variant', 'splice_donor_variant&coding_sequence_variant', 'stop_lost&NMD_transcript_variant', 'frameshift_variant&splice_region_variant', 'frameshift_variant&splice_region_variant&NMD_transcript_variant', 'stop_gained&NMD_transcript_variant', 'stop_retained_variant', 'inframe_deletion&NMD_transcript_variant', 'stop_gained&splice_region_variant', 'frameshift_variant&NMD_transcript_variant', 'stop_lost', 'start_retained_variant&5_prime_UTR_variant', 'coding_sequence_variant', 'frameshift_variant&stop_lost', 'frameshift_variant&stop_retained_variant', 'protein_altering_variant', 'start_lost&splice_region_variant', 'start_lost&splice_region_variant&NMD_transcript_variant', 'frameshift_variant&splice_region_variant&intron_variant', 'stop_lost&splice_region_variant&NMD_transcript_variant', 'stop_retained_variant&NMD_transcript_variant', 'splice_donor_variant&coding_sequence_variant&intron_variant', 'splice_donor_variant&coding_sequence_variant&intron_variant&NMD_transcript_variant', 'inframe_deletion&splice_region_variant', 'stop_gained&frameshift_variant&NMD_transcript_variant', 'stop_gained&frameshift_variant', 'stop_lost&splice_region_variant', 'splice_acceptor_variant&coding_sequence_variant', 'splice_acceptor_variant&coding_sequence_variant&NMD_transcript_variant', 'inframe_insertion&NMD_transcript_variant', 'stop_gained&splice_region_variant&NMD_transcript_variant', 'splice_acceptor_variant&coding_sequence_variant&intron_variant&NMD_transcript_variant', 'protein_altering_variant&NMD_transcript_variant', 'stop_retained_variant&3_prime_UTR_variant', 'inframe_insertion&stop_retained_variant', 'incomplete_terminal_codon_variant&coding_sequence_variant', 'frameshift_variant&stop_retained_variant&NMD_transcript_variant', 'frameshift_variant&start_lost', 'splice_donor_variant&frameshift_variant', 'coding_sequence_variant&3_prime_UTR_variant', 'inframe_insertion&splice_region_variant', 'inframe_insertion&splice_region_variant&NMD_transcript_variant', 'splice_acceptor_variant&frameshift_variant', 'coding_sequence_variant&NMD_transcript_variant', 'frameshift_variant&start_lost&start_retained_variant', 'splice_acceptor_variant&coding_sequence_variant&3_prime_UTR_variant&intron_variant', 'splice_acceptor_variant&coding_sequence_variant&3_prime_UTR_variant&intron_variant&NMD_transcript_variant', 'start_lost&5_prime_UTR_variant', 'stop_gained&frameshift_variant&splice_region_variant', 'stop_lost&3_prime_UTR_variant']

liste_codant_SNP = ['missense_variant', 'splice_region_variant', 'frameshift_variant', 'start_lost', 'inframe_deletion', 'splice_acceptor_variant', 'stop_gained', 'inframe_insertion', 'splice_donor_variant', 'stop_lost', 'protein_altering_variant',"initiator_codon_variant"]


##################################
#list of consequence of mismatchs#
##################################

counted_cons = set(["start_lost","start_retained_variant","protein_altering_variant", "stop_gained", "stop_lost", "initiator_codon_variant", "inframe_insertion", "inframe_deletion", "missense_variant", "frameshift_variant","splice_donor_variant","splice_acceptor_variant","splice_region_variant"])
#counted_cons = set(["start_lost","start_retained_variant","stop_gained", "stop_lost", "initiator_codon_variant", "inframe_insertion", "inframe_deletion", "missense_variant", "frameshift_variant"])
#counted_cons = set(["missense_variant"])

splice_counted_cons=set(["splice_donor_variant","splice_acceptor_variant","intron_variant","splice_region_variant"])


splice_only_counted_cons=["splice_donor_variant","splice_acceptor_variant","splice_region_variant"]

##################################################
#list of consequences of mismtachs no interesting#
##################################################

no_counted_cons = set(["stop_retained_variant"])

def select_with_ratio_median(nt_AD, type_nt,read,median):
    ''' documentation in preparation '''
    
    if len(type_nt)!=len(read):

        read.remove(min(read))

    
    ratio=[]
    majoritaire=0

    DP=sum(read)

    for i in read:
        
        ratio.append(i/DP)
        
        
    for i in range(len(ratio)):
        
        if ratio[i]>0.2+median:
            
            majoritaire=i

    if majoritaire>0:

        nt_AD=str(nt_AD[i])
        read=str(read[i])
        type_nt=str(type_nt[i])

    return nt_AD,type_nt,read




def file_ensembl(path):
    ''' documentation in preparation '''
   
    dic_retour={}
    dic_entete={}
    
    list_file=os.listdir(path)
    
    for i in range(len(list_file)):
        
        files=open(path+list_file[i])
        
        read=files.readline()
        
        while read:
            
            row=read[:-1].split(",")
            
            if "Ensembl Gene ID" in row[0] or "Gene stable ID" in row[0]:

                for j in range(2,len(row)):

                    dic_entete[j]=row[j]
                    
            else:

                gene=row[0]
                transcrit=row[1]
                
                if gene not in dic_retour:

                    dic_retour[gene]={}
                    
                if transcrit not in dic_retour[gene]:    

                    dic_retour[gene][transcrit]=[]
                
                for j in dic_entete:

                    if row[j] !="":
                    
                        dic_retour[gene][transcrit].append(row[j])      
            
            read=files.readline()
        
        
        files.close()
        dic_entete={}
        
    return dic_retour

def extract_best_cons(liste_cons):
    ''' documentation in preparation '''
    liste_cons_2=[]
    
    
    for info in liste_cons:
        
        cons=info.split("|")[3]
        
        liste_cons_2.append(cons)

    
    count_cons=0
    best_cons=""
    
    for cons in liste_cons_2:

        if liste_cons_2.count(cons) > count_cons:
            
            best_cons=cons
            count_cons=liste_cons_2.count(cons)
    
    return best_cons

def extract_best_gene(liste_cons):
    ''' documentation in preparation '''
    liste_cons_2=[]
    
    
    for info in liste_cons:
        
        cons=info.split("|")[1]
        
        liste_cons_2.append(cons)

    count_cons=0
    best_cons=""
    
    for cons in liste_cons_2:

        if liste_cons_2.count(cons) > count_cons:
            
            best_cons=cons
            count_cons=liste_cons_2.count(cons)
        
    return best_cons




def split_and_filtre_CSQ(info):
    ''' documentation in preparation '''
    CSQ = []
    
    if isinstance(info,list) == False:
        
        info=[info]

    for gene in info:

        consequence = gene.split("|")[3].split("&")

        #there is at least one consequence of interest
        if len(set(counted_cons)&set(consequence)) > 0:

            #there is not consequence of no interest
            if len(set(no_counted_cons)&set(consequence)) == 0:

                CSQ.append(gene)


    return CSQ




    for gene in info:

        consequence = gene.split("|")[3].split("&")

        #there is at least one coding consequence
        if len(set(liste_codant_SNP)&set(consequence)) > 0:

            CSQ.append(gene)

    return CSQ


def choix_nt_AD_AMS(ad, ref, alt):
    #classic choice with an alt and a ref according the min and max AD
    ''' documentation in preparation '''

    retour = []

    type_nt = []
    
    read=[]

    for i in range(len(ad)):
        
        if i == 0:

            retour.append(ref)
            type_nt.append("ref")
            

        else:

            retour.append(alt[i-1])
            type_nt.append("alt")
            
        read.append(ad[i])


    return retour, type_nt,read



    retour = []

    ##########################
    ##standard decision of PL#
    ##########################

    #return the smaller PL and his position
    min_, pos_PL = min_PL(pl)

    #keep the standard decision
    if (min_) <= crit_PL:

        #homozygous ref
        if pos_PL == 0:

            retour.append(ref)

        #heterozygous ref alt
        elif pos_PL == 1:

            retour.append(ref)
            retour.append(alt[0])

        #homozygous alt
        elif pos_PL == 2:

            retour.append(alt[0])

    #not keep the standard decision
    else:

        retour = []

    return retour

def choix_nt_AD_AMS_n(ad, gt, ref, alt,read):

    #classic choice with several alt and one ref
    #according the min and max AD
    ''' documentation in preparation '''
    retour = []
    type_nt = []
    read_retour= []

    max_GT = int(max(gt.split("/")))+1

    for i in range(max_GT):

        if str(i) in gt:

            #add reference
            if i == 0:
                
                retour.append(ref)
                type_nt.append("ref")


            #add alt i
            else:

                retour.append(alt[i-1])
                type_nt.append("alt_"+str(i))

    if len(read)==3:
        
        for i in range(len(type_nt)):
            
            type_nnt=type_nt[i]
            
            if type_nnt=="ref":
                
                read_retour.append(read[0])
                
            else:
                
                read_retour.append(read[i+1])

    else:
        
        read_retour=read


    return retour, type_nt,read_retour



    list_borne = liste_alt_f(len(alt))

    min_, pos_min = min_PL(plo)

    ###########################

    #keep the decision
    if (min_) <= crit_PL:

        for i in range(len(list_borne)):

            min_b = list_borne[i][0]
            max_b = list_borne[i][1]

            #heterozygous with alt and ref
            if pos_min == min_b:

                retour.append(ref)
                retour.append(alt[i])

            #alt is homozygous
            elif pos_min == max_b:

                retour.append(alt[i])

            #two alt heterozygous
            elif pos_min > min_b and pos_min < max_b:

                if len(alt) > 1:

                    #ajout du permier alt
                    retour.append(alt[i])
                    alt_2 = [j for j in alt if j != alt[i]]

                    #ajout de 2 alt
                    pos_alt_2 = pos_min-min_b-1
                    retour.append(alt_2[pos_alt_2])

    #not keep the decision
    else:

        retour = []

    return retour


    retour = []
    
    if isinstance(csqr,list) == False:
        
        csqr=[csqr]

    for sample in csqr:

        gene = sample.split("|")[4]

        if gene not in retour and gene != "":

            retour.append(gene)

    return retour

def detect_synonyme(cons):
    ''' documentation in preparation '''
    retour=False

    for c in cons:
        
        if "synonymous_variant" in c:

            retour=True


    return retour

PARSER = argparse.ArgumentParser()

CHOICE = PARSER.add_argument_group('Choix cohorte')
CHOICE.add_argument('-v', '-path_in', type=str, help='Chemin d acces du fichier compile donneur-receveur', required=True)

CHOICE.add_argument('-d', '-donneur', type=str, help='nom du donneur', required=True)
CHOICE.add_argument('-r', '-receveur', type=str, help='nom du receveur', required=True)
CHOICE.add_argument('-o', '-out_dir', type=str, help='Fichier de sortie', required=False, default='out')

CHOICE.add_argument('-s_d', '-sexe_donneur', type=str, choices=["M", "F"], help='sexe du donneur', required=True)
CHOICE.add_argument('-s_r', '-sexe_receveur', type=str, choices=["M", "F"], help='sexe du receveur', required=True)

CHOICE.add_argument('-CSQ', '-fields', type=str, choices=["CSQ", "VariantEffectPrediction"], help='sexe du receveur', required=True)
CHOICE.add_argument('-pe', '-path_ensembl', type=str, help='Chemin d acces du fichier compile donneur-receveur', required=True)


CHOICE.add_argument('-min', '-min_depth', type=int, help='sexe du receveur', required=True)
CHOICE.add_argument('-max', '-max_depth', type=int, help='Chemin d acces du fichier compile donneur-receveur', required=True)
CHOICE.add_argument('-PL', '-phred', type=float, help='sexe du receveur', required=False)
CHOICE.add_argument('-z', '-zigotie', type=str, help='Chemin d acces du fichier compile donneur-receveur', required=False,default=0.95)

CHOICE.add_argument('-c', '-cr', type=str, help='compensate the reference', required=False,default=False)



CHOICE.add_argument('-rd', '-ratio_doneur', type=str , required=False)
CHOICE.add_argument('-rr', '-ratio_receveur', type=str, required=False)

ARGS = PARSER.parse_args()




compense_reference=ARGS.c

if ARGS.rd:
    
    ratio_donneur_median=float(ARGS.rd)
    
    
    
if ARGS.rr:
    
    ratio_receveur_median=float(ARGS.rr)



if ARGS.v:
    
    
    
    LISTE_CHROMOSOME = []
    
    min_depth = ARGS.min
    max_depth = ARGS.max
    crit_PL = ARGS.PL
    crit_ratio=ARGS.z

    

    path_gene_name=ARGS.pe+"gene_name/"
    
    dic_gene_name=file_ensembl(path_gene_name)
    
    info=ARGS.CSQ

    ######################
    #Retrieving variables#
    ######################

    SEXE_D = ARGS.s_d
    SEXE_R = ARGS.s_r

    NAME_D = ARGS.d
    NAME_R = ARGS.r

    MERGED_NAME = NAME_D+"_"+NAME_R

    PATH_DIR_OUT = ARGS.o

    #################################
    ###info quality and read deepth##
    #################################

    DP_DONNEUR_ALL = []
    DP_RECEVEUR_ALL = []

    DP_DONNEUR_CODANT = []
    DP_RECEVEUR_CODANT = []

    DP_DONNEUR_MISMATCH = []
    DP_RECEVEUR_MISMATCH = []

    #############
    #output file#
    #############

    PATH_OUT = PATH_DIR_OUT+MERGED_NAME+"/temp/info_nt.txt"
    F_OUT = open(PATH_OUT, "w")
    F_OUT.write("##CHROMOSOME\tN_SNP\tN_gene\tInfo\tSNP_Donneur\tSNP_Receveur\tread_donneur\tread_receveur\tQual\ttype_nt_d\ttype_nt_r\tref\talt\tinfo_init\n")


    VCF_READER = vcf.Reader(open(ARGS.v, 'rb'))

    #########################################
    #variables of ratio of read of HLA locus#
    #########################################

    RATIO_D_MAX = []
    RATIO_R_MAX = []

    RATIO_D_MIN = []
    RATIO_R_MIN = []
    
    dic_DP_DONNEUR_ALL_v3 = {}
    dic_DP_RECEVEUR_ALL_v3 = {}
    
    
    ######################
    #Variables for splice#
    ######################
    
    save=0
    liste_gene_save=[]
    

    for record in VCF_READER:
        
        #filtering consequences
        
        CSQ = split_and_filtre_CSQ(record.INFO[info])
        chromosome = str(record.CHROM)
        SNP = str(record.POS)
        
        
        synonyme=detect_synonyme(record.INFO[info])
        
   
        if chromosome not in dic_DP_DONNEUR_ALL_v3:

            dic_DP_DONNEUR_ALL_v3[chromosome]={}
            
            
        if SNP not in dic_DP_DONNEUR_ALL_v3[chromosome]:
                     
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]={}
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]=""
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]=""
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["consequence"]=""
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Codant"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Mismatch"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["DP"]=0
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["synonyme"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_ref"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_alt"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMDR"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMRD"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["GT"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Compense"]=0

  
        if chromosome not in dic_DP_RECEVEUR_ALL_v3:
                
                dic_DP_RECEVEUR_ALL_v3[chromosome]={}
                
                
        if SNP not in dic_DP_RECEVEUR_ALL_v3[chromosome]:
                
            
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]={}
            
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]=""
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]=""
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["consequence"]=""
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Codant"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Mismatch"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["DP"]=0
        
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["synonyme"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_ref"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_alt"]=0 
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMDR"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMRD"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["GT"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Compense"]=0
            



        #to see the script progress
        if chromosome not in LISTE_CHROMOSOME:

            print(chromosome)
            LISTE_CHROMOSOME.append(chromosome)

        ref = "".join(list(record.REF))
        alt = list(map(str,record.ALT))

        
        if "-" in ref:
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_ref"]=1
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_ref"]=1
            
            
        if "-" in alt:
            
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_alt"]=1
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_alt"]=1
        
        
        #synonyme
        
        
        if synonyme:
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["synonyme"]=1
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["synonyme"]=1 
        

        #rs
        rs = record.ID
        
        if record.QUAL:
        
            qual = int(round(record.QUAL, 0))
                
        else:
            
            qual="NA"

        #optional
        filter_ = record.FILTER

        donneur = record.genotype(NAME_D)
        receveur = record.genotype(NAME_R)
        
        dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["GT"]=str(".".join(donneur["GT"].split("/")))
        dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["GT"]=str(".".join(receveur["GT"].split("/")))
        
        read_d = donneur["AD"]
        read_r = receveur["AD"]
        

        #########################################
        #Correction du donneur DP et receveur DP#
        #########################################
        
        
        if donneur["AD"]!=None:

            DP_donneur=sum(donneur["AD"])
            
        else:
            
            DP_donneur=0
            
            
        if receveur["AD"]!=None:

            DP_receveur=sum(receveur["AD"])
            
        else:
            
            DP_receveur=0


        ############################
        #####Selection of reads#####
        ############################

        #######
        #Donor#
        #######

        if donneur["AD"]:
            
            

            #classic choice with one alt
            if len(donneur["AD"]) == 2 and len(alt) == 1:

                nt_donneur_AD, type_nt_d,read_d= choix_nt_AD_AMS(donneur["AD"], ref, alt)
                

            else:

                nt_donneur_AD, type_nt_d,read_d = choix_nt_AD_AMS_n(donneur["AD"], donneur["GT"], ref, alt,read_d)

                


        else:
            
            if compense_reference:
            
            
                nt_donneur_AD=[ref,alt[0]]
                type_nt_d=["ref","alt"]
                read_d=[10,0]
                dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Compense"]=1
                DP_donneur=sum(donneur["AD"])
                
            else:
                
                nt_donneur_AD = []
                type_nt_d = []
                read_d=0
                compense_donneur=0
            



        if donneur["AD"]!=None:
            
            if read_d[0]+read_d[1]>10:
                
                if donneur["GT"].split("/")==['0','0']:
                
                    print(donneur["GT"].split("/"))
                    print(read_d,read_d[0])
                    print("\n")

            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]=read_d[0]
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]=read_d[1]
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["DP"]=dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]+dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]
            
            
        else:
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]=0
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["DP"]=0
            
            

            
            
        best_cons=extract_best_cons(record.INFO[info])
        
        best_gene=extract_best_gene(record.INFO[info])
        

        if best_gene!="":

            
            gene_name=dic_gene_name[best_gene][list(dic_gene_name[best_gene].keys())[0]][0]
            
        else:
            
            gene_name="."
            
 
            
        dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["consequence"]=best_cons
        dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["consequence"]=best_cons
        
        dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["gene"]=gene_name
        dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["gene"]=gene_name   
        

        ##########
        #Receveur#
        ##########

        if receveur["AD"]:
            

            #classic choice with one alt
            if len(receveur["AD"]) == 2 and len(alt) == 1:

                nt_receveur_AD, type_nt_r,read_r = choix_nt_AD_AMS(receveur["AD"], ref, alt)


            else:


                nt_receveur_AD, type_nt_r,read_r = choix_nt_AD_AMS_n(receveur["AD"], receveur["GT"], ref, alt,read_r)


            #######################################
            #Recording of reads for plots coverage#
            #######################################


        #without informations, the read is a ref
        else:

            if compense_reference:
            
            
                nt_receveur_AD=[ref,alt[0]]
                type_nt_r=["ref","alt"]
                read_r=[10,0]
                dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Compense"]=1
                DP_receveur=sum(receveur["AD"])
                
            else:
                
                nt_receveur_AD = []
                type_nt_r = []
                read_r=0
                compense_donneur=0
            
        if receveur["AD"]!=None:


            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]=read_r[0]
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]=read_r[1]
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["DP"]=dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]+dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]
            
        else:
            
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]=0
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]=0   
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["DP"]=0

        #############
        #Compute AMS#
        #############

        if CSQ:
            
            
            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Codant"]=1
            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Codant"]=1

            #############
            #Compute AMS#
            #############
            

            if (len(nt_donneur_AD)!=0 and len(nt_receveur_AD)!=0):

                
                if (DP_donneur >= min_depth and DP_donneur <= max_depth) and ( DP_receveur >= min_depth and DP_receveur <= max_depth): 
                    
                    #############
                    #Compute AMS#
                    #############
                    
                    
                    nt_donneur_AD_init, type_nt_d_init,read_d_init= nt_donneur_AD, type_nt_d,read_d
                    nt_receveur_AD_init, type_nt_r_init,read_r_init=nt_receveur_AD, type_nt_r,read_r
                    

                    
                    
                    if ratio_donneur_median:
                        

                        nt_donneur_AD, type_nt_d,read_d=select_with_ratio_median(nt_donneur_AD, type_nt_d,read_d,ratio_donneur_median)
                        
                    
                    if ratio_receveur_median:
                    
                        nt_receveur_AD, type_nt_r,read_r=select_with_ratio_median(nt_receveur_AD, type_nt_r,read_r,ratio_receveur_median)
                    

                    comp="nt_donneur_AD_init="+",".join(nt_donneur_AD_init)+"|type_nt_d_init="+",".join(type_nt_d_init)+"|read_d_init="+",".join(map(str,read_d_init))+"|nt_receveur_AD_init="+",".join(nt_receveur_AD_init)+"|type_nt_r_init="+",".join(type_nt_r_init)+"|read_r_init="+",".join(map(str,read_r_init))
   
                    mismatch = list(set(nt_donneur_AD) ^ set(nt_receveur_AD))
                    
                    MMDR = list(set(nt_donneur_AD).difference(set(nt_receveur_AD)))
                    MMRD = list(set(nt_receveur_AD).difference(set(nt_donneur_AD)))

                    if mismatch != []:

                        if len(MMDR)>0:
     
                            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMDR"]=1
                            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMDR"]=1
                            
                            
                        if len(MMRD)>0:
    
                            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMRD"]=1
                            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMRD"]=1
                            
                            
                            dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Mismatch"]=1
                            dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Mismatch"]=1 
                        
    
                        ############################
                        #Data preparation for write#
                        ############################
    
                        if type(read_d) == list:
    
                            read_d = ",".join(map(str, read_d))
    
                        if not read_d:
    
                            read_d = "None"
    
    
                        if type(read_r) == list:
    
                            read_r = ",".join(map(str, read_r))
    
                        if not read_r:
    
                            read_r = "None"
    
                        ########################################
                        #adjustement according the sex patients#
                        ########################################
    
                        if chromosome == "chrY":
    
                            if SEXE_D == "F":
    
                                read_d = "None_sex"
                                type_nt_d = []
    
                            if SEXE_R == "F":
    
                                read_r = "None_sex"
                                type_nt_r = []
    
                        if read_d != "None_sex" and read_r != "None_sex":
                            

                            F_OUT.write(chromosome+"\t"+SNP+"\t"+"gene"+"\t"+",".join(CSQ)+"\t"+",".join(nt_donneur_AD)+"\t"+",".join(nt_receveur_AD)+"\t"+read_d+"\t"+read_r+"\t"+str(qual)+"\t"+",".join(type_nt_d)+"\t"+",".join(type_nt_r)+"\t"+ref+"\t"+",".join(alt)+"\t"+comp+"\n")

             
            
    F_OUT.close()
    

    with open(PATH_DIR_OUT+MERGED_NAME+"/temp/donneur_DP_v3.csv" ,"w") as f:
        WRITER = csv.writer(f)
        
        WRITER.writerow(["Chromosome","SNP","Read_1","Read_2","DP","Consequence","Coding","Mismatch","Gene_name","Synonyme","Gap_ref","Gap_alt","ratio_ref","ratio_alt","MMDR","MMRD","GT"])
        
        
        for chromosome in dic_DP_DONNEUR_ALL_v3:
            
            for SNP in dic_DP_DONNEUR_ALL_v3[chromosome]:
                
                if dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]=="":
                    
                    dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]=0
                    
                    
                if dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]=="":
                    
                    dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]=0
                    
                    
                    
                if dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"]+dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]>0:

                    ratio_ref=str(int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"])/(int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"])+int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"])))
                    ratio_alt=str(int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"])/(int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"])+int(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"])))

                else:
                    
                    ratio_ref=0
                    ratio_alt=0
                    
                    
                if dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["DP"]>=10:
                    
    
                    buf=[]
                    
                    buf=[chromosome,SNP,dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R1"],str(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["R2"]),dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["DP"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["consequence"],str(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Codant"]),str(dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Mismatch"]),dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["gene"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["synonyme"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_ref"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["Gap_alt"],ratio_ref,ratio_alt,dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMDR"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["MMRD"],dic_DP_DONNEUR_ALL_v3[chromosome][SNP]["GT"]]
    
                    WRITER.writerow(buf)


    with open(PATH_DIR_OUT+MERGED_NAME+"/temp/receveur_DP_v3.csv" ,"w") as f:

        WRITER = csv.writer(f)
        
        WRITER.writerow(["Chromosome","SNP","Read_1","Read_2","DP","Consequence","Coding","Mismatch","Gene_name","Synonyme","Gap_ref","Gap_alt","ratio_ref","ratio_alt","MMDR","MMRD","GT"])
        
        
        for chromosome in dic_DP_RECEVEUR_ALL_v3:
            
            for SNP in dic_DP_RECEVEUR_ALL_v3[chromosome]:
                
                
                if dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]=="":
                    
                    dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]=0
                    
                    
                if dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]=="":
                    
                    dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]=0
                    
                    
                if dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"]+dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]>0:

                    ratio_ref=str(int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"])/(int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"])+int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"])))
                    ratio_alt=str(int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"])/(int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"])+int(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"])))

                else:
                    
                    ratio_ref=0
                    ratio_alt=0
                    
                if dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["DP"]>=10:

                    buf=[]
                    buf=[chromosome,SNP,dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R1"],str(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["R2"]),dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["DP"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["consequence"],str(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Codant"]),str(dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Mismatch"]),dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["gene"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["synonyme"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_ref"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["Gap_alt"],ratio_ref,ratio_alt,dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMDR"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["MMRD"],dic_DP_RECEVEUR_ALL_v3[chromosome][SNP]["GT"]]
     
                    WRITER.writerow(buf)

  


