#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 17:00:02 2018

@author: clerc
"""


import os
import argparse
import csv
import bisect

from fonctions_constructions import construction_peptide_v2


splice_only_counted_cons=["splice_donor_variant","splice_acceptor_variant","splice_region_variant"]


def filtre_ratio(nt_patient_init,read_patient,type_nt):
    ''' documentation in preparation '''
    nt_patient_init_retour=[]
    read_patient_retour=[]
    type_nt_retour=[]
    
    ratio_all=[]
    tag_ratio=[]

    DP_total=sum(read_patient)
    
    for ratio in read_patient:
        
        ratio_all.append(ratio/DP_total)
        
        if ratio/DP_total >= 0.95:
            
            tag_ratio.append(1)
            
        else :
            
            tag_ratio.append(0)
        
        
    if 1 in tag_ratio:
        
        pos_max=tag_ratio.index(1)
        
        nt_patient_init_retour=[nt_patient_init[pos_max]]
        read_patient_retour=[read_patient[pos_max]]
        
        
        if pos_max==0:
            
            type_nt_retour=["Homozygote_Ref"]
            
        elif pos_max==1:
            
            
            type_nt_retour=["Homozygote_Alt"]
        
        
        
    else:
        
        nt_patient_init_retour=nt_patient_init
        read_patient_retour=read_patient
        type_nt_retour=type_nt  

    return nt_patient_init_retour,read_patient_retour,type_nt_retour

def open_bed_file(path_in):
    ''' documentation in preparation '''    

    dic_bed={}
    dic_bed_v2={}
    
    
    file_in=open(path_in,"r")
    
    read=file_in.readline()
    
    
    while read:
        
        
        row=read[:-1].split("\t")

        chromosome=row[0].replace("chr","")
        SNP_begin=int(row[1])
        SNP_end=int(row[2])
        #gene_name=row[3]

        
        if chromosome not in dic_bed:
            
            dic_bed[chromosome]=[]
            print(chromosome)
            
        if (SNP_begin,SNP_end) not in dic_bed[chromosome]:
            
            dic_bed[chromosome].append(SNP_begin)
            dic_bed[chromosome].append(SNP_end)
            
            
        if chromosome not in dic_bed_v2:
            
            dic_bed_v2[chromosome]={}
            
            
        if (SNP_begin,SNP_end) not in dic_bed_v2[chromosome]:
            
            dic_bed_v2[chromosome][(SNP_begin,SNP_end)]={}

        read=file_in.readline()


    file_in.close()
    
    
    return dic_bed,dic_bed_v2


def choose_type_nt(type_nt_donneur_init,type_nt_receveur_init):
    ''' documentation in preparation '''    
    
    retour=[]
    retour_2=[]

    if len(type_nt_donneur_init) == 1 and len(type_nt_receveur_init)==1:
        
        retour=type_nt_donneur_init
        retour_2=type_nt_receveur_init
        
    elif len(type_nt_donneur_init) == 1 and len(type_nt_receveur_init)==2:
        
        retour = type_nt_donneur_init
        retour_2 = list(set(type_nt_receveur_init) - set(list(set(type_nt_donneur_init) & set(type_nt_receveur_init))))
        

    elif len(type_nt_donneur_init) == 2 and len(type_nt_receveur_init)==1:
 
        retour = list(set(type_nt_donneur_init) - set(list(set(type_nt_receveur_init) & set(type_nt_donneur_init))))
        retour_2 = type_nt_receveur_init

    elif len(type_nt_donneur_init) == 2 and len(type_nt_receveur_init)==2: 

        retour = list(set(type_nt_donneur_init) - set(list(set(type_nt_donneur_init) & set(type_nt_receveur_init))))  
        retour_2 = list(set(type_nt_receveur_init) - set(list(set(type_nt_donneur_init) & set(type_nt_receveur_init))))
    

    return retour,retour_2

def choose_nt(nt_donneur_init,nt_receveur_init):
    ''' documentation in preparation '''    
    nt_donneur=[]
    nt_receveur=[]
    
    if len(nt_donneur_init) == 1 and len(nt_receveur_init)==1:

        nt_donneur=nt_donneur_init
        nt_receveur=nt_receveur_init

        
    elif len(nt_donneur_init) == 1 and len(nt_receveur_init)==2:


        nt_donneur = nt_donneur_init
        nt_receveur = list(set(nt_receveur_init) - set(list(set(nt_donneur_init) & set(nt_receveur_init))))

    elif len(nt_donneur_init) == 2 and len(nt_receveur_init)==1:

    
        nt_donneur = list(set(nt_donneur_init) - set(list(set(nt_receveur_init) & set(nt_donneur_init))))   
        nt_receveur = nt_receveur_init

        
    elif len(nt_donneur_init) == 2 and len(nt_receveur_init)==2:

        
        nt_donneur = list(set(nt_donneur_init) - set(list(set(nt_donneur_init) & set(nt_receveur_init))))  
        nt_receveur = list(set(nt_receveur_init) - set(list(set(nt_donneur_init) & set(nt_receveur_init))))
        
        
    #a finir pour gerer les longueures superieurs a 2
    
    else:
        print("fonction choose nt construct_peptides_3 lg sup a 2 a finir")
        print(nt_donneur_init,nt_receveur_init)
        print("\n")

        
    return nt_donneur,nt_receveur
    

    

def check_gap(nt_donneur,nt_receveur,codon):
    ''' documentation in preparation '''    
    retour=0
    
    
    if "-" in codon:
        
        retour=1
        
    for nt in nt_donneur:
        
        if "-" in nt:
            
            retour=1
            
            
    for nt in nt_receveur:
        
        if "-" in nt:
            
            retour=1   

    return retour
    


def compute_ratio_v2(read):
    ''' documentation in preparation '''    
    
    somme=sum(read)
    
    
    return [read[0]/somme,read[1]/somme]

def compute_ratio(read):
    ''' documentation in preparation '''
    if read=="None":
        
        return read
        
    elif len(read.split(","))==1:
        
        return "None"
        
    elif int(read.split(",")[0])==0 and int(read.split(",")[1])==0:
        
        return "None"
    
    else:

        return str(round(float(read.split(",")[0])/(int(read.split(",")[0])+float(read.split(",")[1])),2))


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
            
                #Construction de l'entete

                for j in range(2,len(row)):

                    dic_entete[j]=row[j]
                    
            else:

                gene=row[0]
                transcrit=row[1]
                
                if (gene in dic_retour)==False:
                    
                    dic_retour[gene]={}
                    
                if (transcrit in dic_retour[gene])==False:
                    
                    dic_retour[gene][transcrit]=[]

                
                for j in dic_entete:

                    if row[j] !="":
                    
                        dic_retour[gene][transcrit].append(row[j])
                    
            
            read=files.readline()
        
        
        files.close()
        dic_entete={}
        
    return dic_retour

def open_exac_table_simplified(path):
    ''' documentation in preparation '''    
    dic_retour={}   
        
    with open(path, 'r') as csvfile:
        
        spamreader = csv.reader(csvfile, delimiter=',')
        
        for row in spamreader:
            
            chromosome=row[0]
            pos=row[1]
            ID=row[2]
            SAS=row[3]
            AMR=row[4]
            AFR=row[5]
            OTH=row[6]
            NFE=row[7]
            FIN=row[8]
            EAS=row[9]
            All=row[10]
            
            if "#" not in chromosome:
            
                if (chromosome in dic_retour)==False:
                        
                        dic_retour[chromosome]={}
    
                if (pos in dic_retour[chromosome])==False:
                         
                    dic_retour[chromosome][pos]={}      
                     
                    dic_retour[chromosome][pos]["ID"]=ID
                    dic_retour[chromosome][pos]["Hz"]={}
             
                    dic_retour[chromosome][pos]["Hz"]["SAS"]=0
                    dic_retour[chromosome][pos]["Hz"]["AMR"]=0
                    dic_retour[chromosome][pos]["Hz"]["AFR"]=0
                    dic_retour[chromosome][pos]["Hz"]["OTH"]=0
                    dic_retour[chromosome][pos]["Hz"]["NFE"]=0
                    dic_retour[chromosome][pos]["Hz"]["FIN"]=0
                    dic_retour[chromosome][pos]["Hz"]["EAS"]=0 
                    dic_retour[chromosome][pos]["Hz"]["All_population"]=0
                    
                dic_retour[chromosome][pos]["Hz"]["SAS"]=SAS
                dic_retour[chromosome][pos]["Hz"]["AMR"]=AMR
                dic_retour[chromosome][pos]["Hz"]["AFR"]=AFR
                dic_retour[chromosome][pos]["Hz"]["OTH"]=OTH
                dic_retour[chromosome][pos]["Hz"]["NFE"]=NFE
                dic_retour[chromosome][pos]["Hz"]["FIN"]=FIN
                dic_retour[chromosome][pos]["Hz"]["EAS"]=EAS
                dic_retour[chromosome][pos]["Hz"]["All_population"]=All
            
            
    return dic_retour


def select_stop(gene,transcrit,fasta_donneur_stop):
    ''' documentation in preparation '''    
    #return a stop sequence if it exists
    
    retour=""
    
    if (gene in fasta_donneur_stop)==True:
        
        if (transcrit in fasta_donneur_stop[gene])==True:
            
            retour=fasta_donneur_stop[gene][transcrit]
            
    return retour



def open_fasta(path_fasta):
    ''' documentation in preparation '''    
    dic_fasta={}
    
    fasta_file=open(path_fasta,"r")
    
    read=fasta_file.readline()
    
    
    while read:

        if ">" in read[0]:
            
            line=read[1:-1].split("|")
            gene=line[0]
            transcrit=line[1]
            transcrit=transcrit.replace("\r","").replace("\n","")
            
            if (gene in dic_fasta)==False:
                
                dic_fasta[gene]={}
            
            if (transcrit in dic_fasta[gene])==False:
                
                dic_fasta[gene][transcrit]=""

        else:

            seq=read.replace("\r","").replace("\n","")
            dic_fasta[gene][transcrit]+=seq

        
        
        read=fasta_file.readline()
    
    fasta_file.close()
    
    return dic_fasta


def compute_AMS_d_r(nt_donneur,nt_receveur):
    ''' documentation in preparation '''
    score=len(list(set(nt_donneur).difference(set(nt_receveur))))

    
    return score
        
def compute_AMS_r_d(nt_donneur,nt_receveur):
    ''' documentation in preparation '''
    score=len(list(set(nt_receveur).difference(set(nt_donneur))))

    
    return score

    

def parser_GTEX(path):
    ''' documentation in preparation '''    
    f_in=open(path,"r")
    
    
    read=f_in.readline()
    
    a=0
    
    dic_GTEX={}
    
    header=[]
    
    while read:

        row=read[:-1].split("\t")
        
        if a==2:
        
            for i in range(2,len(row)):
                
                header.append(row[i])
                
        elif a>2:

            gene=row[0].split(".")[0]
            
            
            if gene not in dic_GTEX:
                
                dic_GTEX[gene]={}
                
            for i in range(2,len(row)-2):
                
                dic_GTEX[gene][header[i]]=""
                
                dic_GTEX[gene][header[i]]=row[i]
    
        a+=1
    
        read=f_in.readline()
    
    
    f_in.close()
    
    return dic_GTEX
    

    
################
#Main Principal#
################


parser = argparse.ArgumentParser()

choice = parser.add_argument_group('Choix cohorte')

choice.add_argument('-e','-path_ensembl',type=str,help='Fichier d entree de la base de donnee Ensembl',required=True)
choice.add_argument('-i','-path_in',type=str,help='Fichier d entree du patient',required=True)
choice.add_argument('-l','-liste_lg',type=str,help='liste des longeures de peptides a construire',required=True)

choice.add_argument('-e_d', '-ethnie_donneur', type=str, help='ethnie du donneur', required=True)

choice.add_argument('-e_r', '-ethnie_receveur', type=str, help='ethnie du receveur', required=True)

choice.add_argument('-gtx', '-gtex', type=str, help='chemin d acces a la base de donnee gtex', required=True)

choice.add_argument('-ex', '-exac', type=str, help='chemin d acces a la base de donnee exac', required=True)

choice.add_argument('-p','-path_peptides',type=str,help='Chemin d acces du repertoire des peptides',required=True)

choice.add_argument('-t','-path_table_premhc',type=str,help='Chemin d acces de la table preMHC',required=True)

choice.add_argument('-b','-bed_file',type=str,help='Chemin d acces de la table preMHC',required=False,default="None")

args = parser.parse_args()


path_bed=args.b

if path_bed != "None":


    dic_bed,dic_bed_v2 = open_bed_file(path_bed)
    
else:
    
    dic_bed={}
    dic_bed_v2={}

liste_longueure=list(map(int,args.l.split(",")))

dic_file_out_d={}
dic_file_out_r={}

path_e_expression=args.e+"expression/organism_part/"
path_e_sequence=args.e+"sequence/"

ethnie_d=args.e_d
ethnie_r=args.e_r

PATH_EXAC=args.ex


#Files sequences

fasta_sequence=open_fasta(path_e_sequence+"sequence.fasta")
fasta_sequence_stop=open_fasta(path_e_sequence+"sequence_stop.fasta")

#ensembl organism expression

dic_organism_part_Ensembl=file_ensembl(path_e_expression)


dic_Hz=open_exac_table_simplified(PATH_EXAC)



dic_GTEX = parser_GTEX(args.gtx) 
    


#path for peptide file for pair of patient
path_patient=args.i
path_out=args.p+"/peptides"




liste_patient=os.listdir(path_patient)
    
#list of consequences unmanaged


#preparation of output files
    
if os.path.isdir(path_out)==False:
    os.mkdir(path_out)
        
if os.path.isdir(path_out+"/receveur")==False:
    os.mkdir(path_out+"/receveur")
        
if os.path.isdir(path_out+"/donneur")==False:
    os.mkdir(path_out+"/donneur")
    
if os.path.isdir(path_out+"/donneur/peptides")==False:
    os.mkdir(path_out+"/donneur/peptides")
    
if os.path.isdir(path_out+"/donneur/info")==False:
    os.mkdir(path_out+"/donneur/info")
    
if os.path.isdir(path_out+"/receveur/peptides")==False:
    os.mkdir(path_out+"/receveur/peptides")
    
if os.path.isdir(path_out+"/receveur/info")==False:
    os.mkdir(path_out+"/receveur/info")
    
    
f_upgrade=open(path_out+"/upgrade_construct_peptides.txt","w")
f_gap=open(path_out+"/gap_non_traite.txt","w")
f_stop=open(path_out+"/seq_stop_a_rajouter.txt","w")

dic_file_out_d={}
dic_file_out_r={}

dic_file_out_d_2={}
dic_file_out_r_2={}

buf_header_d=[]
buf_header_r=[]

for f in range(len(liste_longueure)):
    
    dic_file_out_d[liste_longueure[f]]=open(path_out+"/donneur/peptides/"+str(liste_longueure[f])+".fasta","w")
    dic_file_out_r[liste_longueure[f]]=open(path_out+"/receveur/peptides/"+str(liste_longueure[f])+".fasta","w")
    
    dic_file_out_d_2[liste_longueure[f]]=open(path_out+"/donneur/info/"+str(liste_longueure[f])+"_info.fasta","w")
    dic_file_out_r_2[liste_longueure[f]]=open(path_out+"/receveur/info/"+str(liste_longueure[f])+"_info.fasta","w")
    
    buf_header_d.append(str(liste_longueure[f])+"_r")
    buf_header_r.append(str(liste_longueure[f])+"_d")


buf_header_d.extend(buf_header_r)

#read of input file
f_in=open(path_patient+"/temp/info_nt_filter.txt","r")
read=f_in.readline()

cle=0

#allows to define the procedure for the function of building the script functions_buildings

nb_peptides_donneur={}
nb_peptides_receveur={}

nb_peptides_donneur["total"]=0
nb_peptides_receveur["total"]=0


fname = args.t+"/table_pre_MHC.csv"

files = open(fname, "w")

writer = csv.writer(files, delimiter=',')

buf_header=["#chromosome","SNP","gene","gene_name","hyp_I","hyp_II","hyp_III","nombre_peptide_donneur","nombre_peptide_receveur","ExAC_donneur","ExAC_receveur","ExAC_all_population","value_GTEX","consequence","consequence_simplified","lg_peptide","transm","localisation_ensembl_rein","ratio_read_donneur","ratio_read_receveur","qual","PDB","domaine","info","SNP_donneur_1","SNP_donneur_2","SNP_receveur_1","SNP_receveur_2","read_donneur_1","read_donneur_2","read_receveur_1","read_receveur_2","type_donneur_1","type_donneur_2","type_receveur_1","type_receveur_2","bed_value","Imput_d","Imput_r"]
buf_header.extend(buf_header_d)
     
writer.writerow(buf_header)
    
upgrade=0
gap=0

a=0
origin=""

liste_chromosome=[]

seq_stop_upgrade=0

while read:

    if "#" not in read[0]:
    

        row=read[:-1].split("\t")

        info=row[3].split("|")

        chromosome=row[0]
        SNP=row[1]
        
        SNP_donneur=row[4]
        SNP_receveur=row[5]
        
        qual=row[10]

        gene=info[1]
        transcrit=info[2]
        gene_name=info[13].split("=")[1]
        consequence=info[3].split("&")
        transm=info[10].split("=")[1]
        PDB=info[12].split("=")[1]
        domaine=info[11].split("=")[1]
    
        nt_donneur_init=row[4].split(",")
        nt_receveur_init=row[5].split(",")

        read_donneur=list(map(int,row[6].split(",")))
        read_receveur=list(map(int,row[7].split(",")))

        consequence_simplified=consequence[0]

        type_nt_donneur_init=row[8].split(",")
        type_nt_receveur_init=row[9].split(",")
        
        
        ref_vcf=row[11]
        alt_vcf=row[12]
        
        Imput_d=row[13]
        Imput_r=row[14]

        
        ratio_donneur=compute_ratio_v2(read_donneur)
        ratio_receveur=compute_ratio_v2(read_receveur)


        #filtre DP sur liste_nt_donneur et liste_nt_receveur a faire et read_donneur read_receveur

        if len(nt_donneur_init)==2:
        
        
            nt_donneur_init,read_donneur,type_nt_donneur_init=filtre_ratio(nt_donneur_init,read_donneur,type_nt_donneur_init)

        
        if len(nt_receveur_init)==2:
        
        
            nt_receveur_init,read_receveur,type_nt_receveur_init=filtre_ratio(nt_receveur_init,read_receveur,type_nt_receveur_init)

        liste_nt_donneur=nt_donneur_init
        liste_nt_receveur=nt_receveur_init

        
        type_donneur_1,type_receveur_1=choose_type_nt(type_nt_donneur_init,type_nt_receveur_init)
        
    
        
        if len(type_donneur_1)==2:
            
            
            type_donneur_1=type_donneur_1[0]
            type_donneur_2=type_donneur_1[1]
            
        elif len(type_donneur_1)==1 :
            
            type_donneur_1=type_donneur_1[0]
            type_donneur_2="NA"
            
        else :
            
            type_donneur_1="NA"
            type_donneur_2="NA"
            
        
        if len(type_receveur_1)==2:
            
            
            type_receveur_1=type_receveur_1[0] 
            type_receveur_2=type_receveur_1[1]
            
        elif len(type_receveur_1)==1 :
            
            type_receveur_1=type_receveur_1[0]
            type_receveur_2="NA"

        else :
            
            type_receveur_1="NA"
            type_receveur_2="NA"
    
        if len(liste_nt_donneur)==2:
        
            SNP_donneur_1=liste_nt_donneur[0]
            SNP_donneur_2=liste_nt_donneur[1]
            
        else:
            
            SNP_donneur_1=liste_nt_donneur[0]
            SNP_donneur_2="NA"
        
        
        if len(liste_nt_receveur)==2:
        
            SNP_receveur_1=liste_nt_receveur[0]
            SNP_receveur_2=liste_nt_receveur[1]
            
        else:
            
            SNP_receveur_1=liste_nt_receveur[0]
            SNP_receveur_2="NA"   


        if len(read_donneur)==2:

            read_donneur_1=int(read_donneur[0])
            read_donneur_2=int(read_donneur[1])
            
        else:
            
            read_donneur_1=int(read_donneur[0])
            read_donneur_2=0      
            
        if len(read_receveur)==2:          
        
            read_receveur_1=int(read_receveur[0])
            read_receveur_2=int(read_receveur[1])
            
        else:
            
            
            read_receveur_1=int(read_receveur[0])
            read_receveur_2=0       

        
        DP_donneur=read_donneur_1+read_donneur_2
        DP_receveur=read_receveur_1+read_receveur_2
        

        hyp_I=0
        hyp_II=0
        hyp_III=0
        
        origin=""


        #Reduction des SNPS pour ne garder que les mismatchs non en commun entre le donneur et le receveur
          
        nombre_peptide_donneur=0
        nombre_peptide_receveur=0
        
        ############################
        #selection des nt d'interets
        
        

        #Permet de compter les scores selon les hypotheses
        hyp_I=compute_AMS_d_r(nt_donneur_init,nt_receveur_init)
        hyp_II=compute_AMS_d_r(nt_donneur_init,nt_receveur_init)
        hyp_III=compute_AMS_r_d(nt_donneur_init,nt_receveur_init)
        

        
        #mise a niveau cu chromosome
        if "chr" in chromosome:

            
            chromosome=chromosome.replace("chr","")
            
            
        if chromosome not in liste_chromosome:
            
            liste_chromosome.append(chromosome)
        
        
        #extraction des donnees de frequences
        if chromosome in dic_Hz:
            
            if SNP in list(dic_Hz[chromosome].keys()):

                Exac_all=str(dic_Hz[chromosome][SNP]["Hz"]["All_population"])
                   
                if len(liste_nt_donneur)==1 and len(liste_nt_receveur)==1:

                    ref=ref_vcf[0]
                    alt=alt_vcf[0]

                    if liste_nt_donneur[0]==ref and liste_nt_receveur[0]==alt:
                        
                        ExAC_receveur=str(Exac_all)
                        ExAC_donneur=str(1-float(Exac_all))
                        
                    elif liste_nt_donneur[0]==alt and liste_nt_receveur[0]==ref:
                    
                        ExAC_donneur=str(Exac_all)
                        ExAC_receveur=str(1-float(Exac_all))
                        
                        
                    else:

                        ExAC_donneur="NA"
                        ExAC_receveur="NA"
                        
                        Exac_all="NA"

                        
                else:
                    
                    ExAC_donneur="NA"
                    ExAC_receveur="NA"
                    
                    Exac_all="NA"
                
            else:
                           
                ExAC_donneur="NA"
                ExAC_receveur="NA"
                
                Exac_all="NA"
                
        else:
                       
            ExAC_donneur="NA"
            ExAC_receveur="NA"
            
            Exac_all="NA"                                               
                                                
        if gene in dic_GTEX:

            value_GTEX=str(dic_GTEX[gene]["Kidney - Cortex"])
            
        else:
            
            value_GTEX="NA"
            
            

        #1
        #on ne construit pas les splices dans un premier temps
        if len(list(set(consequence)&set(splice_only_counted_cons)))==0:

            #2
            #on ne construit les peptides que pour les sequences genetiques connues
            if (gene in fasta_sequence)==True:    
                        
                #3
                if (transcrit in fasta_sequence[gene])==True:
         
                    #4
                    if fasta_sequence[gene][transcrit]!="Sequence unavailable":
                                    
                        sequence=fasta_sequence[gene][transcrit]
    
                        if "-" not in info[5] and info[5]!="":
    
                            pos_sequence=int(info[5])
                                        
                        else:
                                                                
                            pos_sequence=info[5]

                        chg_codon=info[8]
                        
                 
                        #selection of stop sequence
                        if "stop_lost" in consequence or "frameshift_variant" in consequence:
    
                            seq_stop=select_stop(gene,transcrit,fasta_sequence_stop)
       
                        else:
                                        
                            seq_stop=""

                        gap=0
                        gap=check_gap(nt_donneur_init,nt_receveur_init,info[8])
                        
                        origin=""

                        #5
                        if gap==0:

                            #building peptides of donor
                            for n in range(len(liste_nt_donneur)):
                                
                                nt=liste_nt_donneur[n]

                                #6
                                if "-" not in nt and "-" not in chg_codon and "?" not in str(pos_sequence):
     
                                    for lg in liste_longueure:

                                        peptide,pos_snp,ev,seq_stop_upgrade,even_choix=construction_peptide_v2(str(nt),lg,sequence,consequence,pos_sequence,chg_codon,seq_stop,gene,transcrit,"|".join(info))

                                        if seq_stop_upgrade==1:
                                            
                                            f_stop.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"\n\n")
                                            seq_stop_upgrade=0

                                        if ev=="upgrade":
                   
                                            upgrade=1
                                            origin=origin+"&ev_fct_contruction_donneur&"
           
                                        else:
                        
                                            if str(peptide)!="" and len(peptide)>=lg:
                                                                                 
                                                dic_file_out_d[lg].write(">"+hex(cle)+"\n")
                                                dic_file_out_d[lg].write(str(peptide)+"\n")
            
                                                
                                                dic_file_out_d_2[lg].write(">"+hex(cle)+"\n")
                                                dic_file_out_d_2[lg].write(chromosome+"|"+SNP+"|"+gene+"|"+"|".join(info[2:len(info)]).replace("/",".")+"|position_mismatch_peptidique="+str(pos_snp)+"|lg_sequence="+str(len(sequence))+"|ev="+str(ev)+"|nt_interet="+nt+"|hyp_I_hyp_II="+str(hyp_I)+"|hyp_III="+str(hyp_III)+"|nt_init_donneur="+",".join(nt_donneur_init)+";nt_init_receveur="+",".join(nt_receveur_init)+"\n")
                                                cle+=1

                                 
                                                nombre_peptide_donneur+=(len(peptide)-lg)+1

                                                if nombre_peptide_donneur<0:
                                                    
                                                    nombre_peptide_donneur=0   
                                                
                                                nb_peptides_donneur["total"]+=nombre_peptide_donneur
                                                
                                                if (lg in nb_peptides_donneur)==False:
                                                    
                                                    nb_peptides_donneur[lg]=0
                                                    
                                                nb_peptides_donneur[lg]=(len(peptide)-lg)+1
                                            
        
                                #6
                                else:                                     
                                                
                                    upgrade=1
                                    origin=origin+"-_ou_?_donneur&"
                                    f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix="+even_choix+"\n")

        
        
                            #building peptides of donor      
                            for n in range(len(liste_nt_receveur)):
                                        
                                nt=liste_nt_receveur[n]

                                #7
                                if "-" not in nt and "-" not in chg_codon and "?" not in str(pos_sequence):
                                    
                                    for lg in liste_longueure:
                                        
                                        peptide,pos_snp,ev,seq_stop_upgrade,even_choix=construction_peptide_v2(str(nt),lg,sequence,consequence,pos_sequence,chg_codon,seq_stop,gene,transcrit,"|".join(info)) 
            
            
                                        if seq_stop_upgrade==1:
                                            
                                            f_stop.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"\n\n")
                                            seq_stop_upgrade=0
            
                                        if ev=="upgrade":
                                            
                                            upgrade=1
                                            origin=origin+"&ev_fct_contruction_receveur&"
                                                
                                        else:
                                            
         
                                            if str(peptide)!="" and len(peptide)>=lg:
                                                
            
                                                dic_file_out_r[lg].write(">"+hex(cle)+"\n")
                                                dic_file_out_r[lg].write(str(peptide)+"\n")
                                             
                                                dic_file_out_r_2[lg].write(">"+hex(cle)+"\n")
                                                dic_file_out_r_2[lg].write(chromosome+"|"+SNP+"|"+gene+"|"+"|".join(info[2:len(info)]).replace("/",".")+"|position_mismatch_peptidique="+str(pos_snp)+"|lg_sequence="+str(len(sequence))+"|ev="+str(ev)+"|nt_interet="+nt+"|hyp_I_hyp_II="+str(hyp_I)+"|hyp_III="+str(hyp_III)+"|nt_init_donneur="+",".join(nt_donneur_init)+";nt_init_receveur="+",".join(nt_receveur_init)+"\n")
                                                        
                                                cle+=1
                                                
                                                
                                                nombre_peptide_receveur+=(len(peptide)-lg)+1
                                                
                                                if nombre_peptide_receveur<0:
                                                    
                                                    nombre_peptide_receveur=0   
                                                
                                                nb_peptides_receveur["total"]+=nombre_peptide_receveur
                                                
                                                if (lg in nb_peptides_receveur)==False:
                                                    
                                                    nb_peptides_receveur[lg]=0
                                                    
                                                nb_peptides_receveur[lg] = (len(peptide)-lg)+1
                                                
                                #7               
                                else: 
                                
                                    upgrade=1
                                    origin=origin+"-_ou_?_receveur&"
                                    even_choix="NA"
                                    
                                    f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix="+even_choix+"\n")

                                    upgrade=0   
                                    
                             
                        #5
                        else:

                            f_gap.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"\n")
                            gap=0
                                    
                    else:
                        

                        f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix=seq_non_dispo\n")

                                    
                else:

                    f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix=seq_transcrit_non_dispo\n")
                                    
            else:

                f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix=seq_gene_non_dispo\n")
                         
        else:

            f_upgrade.write(">"+"|".join(info)+"|nt_donneur="+".".join(nt_donneur_init)+"|nt_receveur="+".".join(nt_receveur_init)+"|origin="+origin+"|even_choix=splice_non_gerer\n")                  
            
 
        ####################
        ###Parser Ensembl###
        ####################
        
        if gene in dic_organism_part_Ensembl:
            
            if transcrit in dic_organism_part_Ensembl[gene]:
                
                if dic_organism_part_Ensembl[gene][transcrit]!=[]:
    
                    localisation_ensembl=",".join(dic_organism_part_Ensembl[gene][transcrit])
                    
                else:
                    
                    localisation_ensembl="NA"

            else:
                
                localisation_ensembl="NA"  

        else:
            
            localisation_ensembl="NA"

            
            
        if "kidney" not in localisation_ensembl:

            localisation_ensembl="NA"
            
        else:

            localisation_ensembl="T"
        
        
        ########################
        ###Fin parser Ensembl###
        ########################

        
    
        buf_donneur=[]
        buf_receveur=[]
        
        for lg in liste_longueure:

            if lg in nb_peptides_donneur:
            
                buf_donneur.append(str(nb_peptides_donneur[lg]))
                
            else:
                
                buf_donneur.append(str(0))
                
                
                
            if lg in nb_peptides_receveur:
            
                buf_receveur.append(str(nb_peptides_receveur[lg]))
                
            else:
                
                buf_receveur.append(str(0))
                
            

        
        ##############################
        #Bed file for common interval#
        ##############################
        
        if path_bed != "None":
            
            if chromosome in dic_bed:
        
                test = bisect.bisect_left(dic_bed[chromosome], int(SNP))
                
    
                #Si le test fonctionne, le SNP est dans un intervalle couvert par le bed commun
                if test%2==1 and test!=0: #test!=len(dic_bed[Chromosome])
    
                    bed_value=1
                    
                else:
                
                    bed_value=0
                
            else:
                
                bed_value=0
            
        else:
            
            bed_value=0

        
        temp=[]      
        
        ratio_donneur=ratio_donneur[0]
        ratio_receveur=ratio_receveur[0]
  
        temp=[chromosome,SNP,gene,gene_name,hyp_I,hyp_II,hyp_III,nombre_peptide_donneur,nombre_peptide_receveur,ExAC_donneur,ExAC_receveur,Exac_all,value_GTEX,".".join(consequence),consequence_simplified,(".").join(map(str,liste_longueure)),transm,localisation_ensembl,str(ratio_donneur),str(ratio_receveur),qual,PDB,domaine,"|".join(info),SNP_donneur_1,SNP_donneur_2,SNP_receveur_1,SNP_receveur_2,read_donneur_1,read_donneur_2,read_receveur_1,read_receveur_2,type_donneur_1,type_donneur_2,type_receveur_1,type_receveur_2,bed_value,Imput_d,Imput_r]
        
  
        temp = temp + buf_donneur + buf_receveur
        
        writer.writerow(temp)


        nombre_peptide_donneur=0
        nombre_peptide_receveur=0
        
        buf_donneur=[]
        buf_receveur=[]
        
     
        nb_peptides_donneur={}
        nb_peptides_receveur={}
        
        nb_peptides_donneur["total"]=0
        nb_peptides_receveur["total"]=0

    read=f_in.readline()  
    

#closing files
for l in dic_file_out_d:
            
    dic_file_out_d[l].close()
    dic_file_out_d_2[l].close()
            
for l in dic_file_out_r:
            
    dic_file_out_r[l].close()
    dic_file_out_r_2[l].close()
    
        
dic_file_out_d={}
dic_file_out_r={}
    
f_in.close()
f_upgrade.close()
f_gap.close()
f_stop.close()