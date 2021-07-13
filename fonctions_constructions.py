#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 16:30:24 2018

@author: clerc
"""


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



def lg_SNP_compute(pos_SNP):
    #return the length of SNP 
    ''' documentation in preparation '''
    un,deux=list(map(int,pos_SNP.split("-")))
    
    return deux-un+1

def separator(chg_codon): 
    #separate both bornes 
    ''' documentation in preparation '''
    chg_codon = (chg_codon.replace("/", " "))
    
    return chg_codon.split(" ")

def SNPposition_to_consider(pos_SNP): 
    #return the lower SNP position
    ''' documentation in preparation '''    
    return int(pos_SNP[:pos_SNP.find("-")])

def only_gap(codon):
    #know if a codon have only gap
    ''' documentation in preparation '''    
    retour=0
    
    if codon.count("-") == len(codon):
        
        retour=1

    return retour

def count_maj(char):
    ''' documentation in preparation '''    
    retour=0
    
    for i in char:
        
        if ord(i)>=65 and ord(i)<=90:
            
            retour+=1
    
    return retour

def nb_minuscule_amont(alt):
    """
    return the number of minuscule upstream
    """
    ''' documentation in preparation '''    
    retour=0
    a=0
    
    for char in alt:

        if a==0:
        
            if char.islower():
                
                retour+=1
                
            else:
                
                a=1

    return retour


def enlever_minuscule(Biomart_codon):

    #remove lowercases of a str for keep only uppercase
    ''' documentation in preparation '''    
    retour=""
    
    for i in range(len(Biomart_codon)):
  
        if ord(Biomart_codon[i])>=65 and ord(Biomart_codon[i])<=90 :
            
            retour+=Biomart_codon[i]

    return retour


def cut_seq(seq_stop,pos_next_stop):
    #cut a sequence according to a psition of SNP
    ''' documentation in preparation '''


    return seq_stop[:int(pos_next_stop)*3]

def det_stop(seq_stop):
    """
    det_stop
    """
    #from a sequence, return a position of first stop codon find
    #return 0 if no stop codon is in the sequence
    ''' documentation in preparation '''
    a=0
    retour=0
    seq_stop_upgrade=0
  
    ll=list(range(0,len(seq_stop),3))
 
    for i in ll:
        
        if a==0:

            if Seq(seq_stop[i:i+3], IUPAC.unambiguous_dna).translate()=="*":
                
                a=1

                retour=(i/3)+1

    if retour>6000:
        
        seq_stop_upgrade=1

    return retour,seq_stop_upgrade

def pos_ajust_alt_2(alt,nt):
    ''' documentation in preparation '''   
    """
    from a read and an alt, search if alt is in nt and how many nucleotides in 
    the nt are upstream before the common sequence with alt
    only possible without gap
    """
    
    nt_alt=0
    pos_aj=0
    
    alt_upper=enlever_minuscule(alt)

    if alt_upper in nt:
        
        if len(nt) >= len(alt_upper):

            nt_alt=1

            
    #return the number of lower case upstream in alt
    pos_aj=nb_minuscule_amont(alt)

    return pos_aj,nt_alt

def pos_ajust_alt(alt,nt):
    """
    from a read and an alt, search if alt is in nt and how many nucleotides in 
    the nt are upstream before the common sequence with alt
    only possible without gap and if the ref( or alt) is a gap
    """
    ''' documentation in preparation '''   
    nt_alt=0

    if alt in nt:
        
        nt_alt=1

    return nt_alt

def choix_particulier_missense(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence):
    """
    choix_particulier_missense
    """
    ''' documentation in preparation '''   
    #we consider a read error and apply the event normally
    if len(nt)>1 and "-" not in chg_codon and "-" not in nt and "-" not in str(pos_sequence):
        
        ev+="&missense_error_read"
        even_choix="chgt_simple"
        
        
    elif len(nt)==1 and "-" not in chg_codon and "-" not in nt and "-"  in str(pos_sequence):
        
        ev+="&missense_error_read_2"
        even_choix="chgt_simple"
        
        pos_sequence=int(str(pos_sequence).split("-")[0])
        
                
    elif len(nt)>1 and "-" not in chg_codon and "-" not in nt and "-"  in str(pos_sequence):
        
        ev+="&missense_error_read_4"
        even_choix="chgt_simple"
        
        pos_sequence=int(str(pos_sequence).split("-")[0])
        
    else:    

        even_choix="upgrade_choix_particulier_missense"
        ev="upgrade"
        
    
    return  even_choix,pos_sequence,nt,lg_deletion,ev

def choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info):
    """
    choix_particulier_insertion
    """    
    ''' documentation in preparation '''

    
    if chg_codon!="":
    
    
        if "-"  in str(pos_sequence) and "-" not in nt and only_gap(separator(chg_codon)[0]):
            

            ##########################
            #preparation of variables#
            ##########################
            
            ref,alt=separator(chg_codon)
            
            #compute the lg of event according to the position SNP reference
            lg_even=lg_SNP_compute(pos_sequence)
            
            #first position in the reference codon  123-125 <- 123
            first_pos=SNPposition_to_consider(pos_sequence)
            
            first_decoup_seq=first_pos+lg_even-1
            
            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
            
            nt_alt=pos_ajust_alt(alt,nt)
    
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################  
            
            seq_alt=sequence[first_pos:first_pos+len(alt)]
            
    
            if alt==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0
    
    
                
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
                
                ev+="&nothing_insert"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
                
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1:
                #ok definitif
    
                ev+="&deletion_ajustee"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos
                lg_deletion=len(seq_alt)
    
                
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0:
                #ok definitif
                ev+="&new_insert_spe"
                even_choix="insert"           
                pos_sequence=first_pos
    
                #warning fir insertion, it'possible to find nucleotides in too much in nt compared to alt
                #we need remove nucleotides in too much
                
                if alt in nt:
                    
                    nt=alt
                    
                else:
    
                    even_choix="upgrade_insertion_1"
                    ev="upgrade"
             
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
                            
                ev+="&nothing_insert_empty"
                even_choix="nothing_insert"
                pos_sequence=first_pos
              
        elif "-"  in str(pos_sequence) and "-" not in nt and not only_gap(separator(chg_codon)[0]):
            

    
            ##########################
            #preparation of variables#
            ##########################
            
            ref,alt=separator(chg_codon)
    
            #compute the lg of event in terms of position of SNP
            lg_even=lg_SNP_compute(pos_sequence)
            
            #first position in the reference codon  123-125 <- 123
            first_pos=SNPposition_to_consider(pos_sequence)
        
            first_decoup_seq=first_pos+lg_even-1

            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
    
            #renvoi le nombre de nt min en amont de l'alt et un bool si l'alt est inclus dans le nt, auquel cas le nt est bien une insertion
            pos_ajust,nt_alt=pos_ajust_alt_2(alt,nt)
            
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################  
   
            seq_alt=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(alt)-1-pos_ajust]

            #test si la sequence alt existe deja dans la sequence de reference
            if alt.upper()==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0   

    
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
    
                ev+="&nothing_insert_2"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
             
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1 :
                
                ev+="&deletion_ajustee_2"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos+pos_ajust-1
    
                lg_deletion=count_maj(alt)
                
    
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0 :
               
                ev+="&new_insert_spe_2"
                even_choix="insert"
    
                pos_sequence=first_pos
    
                if enlever_minuscule(alt) in nt:
                    
                    nt=enlever_minuscule(alt)
                    
                else:
                    
                    even_choix="upgrade_insertion_2"
                    ev="upgrade"
            
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
    
    
                ev+="&nothing_insert_empty_2"
                even_choix="nothing_insert"
                pos_sequence=first_pos
                
                
        elif "-" not in str(pos_sequence) and "-" not in nt and "-" not in chg_codon:
            
    
            ###########################
            #preparation des variables#
            ###########################
            
            ref,alt=separator(chg_codon)
            
            #compute the lg of event in terms of position of SNP
            lg_even=len(ref)-len(alt)
            
    
            #first position in the reference codon  123-125 <- 123
            first_pos=pos_sequence
            
            first_decoup_seq=first_pos+lg_even-1
            
            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
    
            pos_ajust,nt_alt=pos_ajust_alt_2(alt,nt)
            
    
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################    
            
            seq_alt=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(alt)-1-pos_ajust]
            
    
            if alt.upper()==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0   
    
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
    
                ev+="&nothing_insert_3"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
             
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1 :
                
                ev+="&deletion_ajustee_3"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos+pos_ajust-1
    
                lg_deletion=count_maj(alt)
                
    
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0 :
               
                ev+="&new_insert_spe_3"
                even_choix="insert"
    
                pos_sequence=first_pos
    
                if enlever_minuscule(alt) in nt:
                    
                    nt=enlever_minuscule(alt)
                    
                else:
                    
                    even_choix="upgrade_insertion_3"
                    ev="upgrade"
            
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
    
    
                ev+="&nothing_insert_empty_3"
                even_choix="nothing_insert"
                pos_sequence=first_pos

        else : 
            
            even_choix="upgrade_insertion_4"
            ev="upgrade"

                
    elif chg_codon=="" and len(nt)==1:
        
        
        if "-" in str(pos_sequence):
        
            lg_even=lg_SNP_compute(pos_sequence)
            first_pos=SNPposition_to_consider(pos_sequence)
            
        else:
            
            lg_even=1
            first_pos=int(pos_sequence)

        nt_courant=sequence[first_pos-1]
        
        nt_suivant=sequence[first_pos]
        
 
        if nt==nt_courant:
            
            ev+="&chg_codon_empty_nothing"
            even_choix="nothing_delet"
            pos_sequence=first_pos
  
        elif nt==nt_suivant:
            
            ev+="&deletion_chg_codon_empty"
            even_choix="deletion_ajustee"
            pos_sequence=first_pos
            lg_deletion=len(nt)
            
        else : 
    
            even_choix="upgrade_insertion_5"
            ev="upgrade"   
        
        
            
    else : 

        even_choix="upgrade_insertion"
        ev="upgrade"   

    return even_choix,pos_sequence,nt,lg_deletion,ev
    
    
def choix_particulier_insertion_v2(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info,type_nt):
    """
    choix_particulier_insertion
    """    
    ''' documentation in preparation '''      
    
    #Preparation des variables
    
    if "-" in str(pos_sequence):
        
        first_pos=int(pos_sequence.split("-")[0])
        second_pos=int(pos_sequence.split("-")[1])
        lg_deletion=second_pos - first_pos
        
    else:
        
        first_pos=int(pos_sequence)
        lg_deletion=0
        second_pos=first_pos
        
        
    lg_nt=len(nt)

    
    if chg_codon != "":
        
        ref=chg_codon.split("/")[0]
        alt=chg_codon.split("/")[1]
        
        ref_upper=ref.upper()
        alt_upper=alt.upper()
        

                
        shift=sequence[first_pos-1:first_pos+lg_deletion]
        
        premier=shift[0]
        dernier=shift[len(shift)-1]
        suivant=sequence[first_pos+lg_deletion]
        
        p=0
        d=0
        s=0
        
        #C'est une insertion et elle est deja faite, donc pas de changemenr
        if nt==shift:
            
            ev+="&nothing_insert"
            even_choix="nothing_insert"
            
        else:

            
            if nt==premier:
                
                p=1
                
            
            if nt==dernier:
    
                d=1
                
                
            if nt==suivant:
                s=1
                
            if len(nt)>1:
                
                if (nt in shift) or (shift in nt):
                    
                    print("a",nt,lg_nt,first_pos,second_pos,lg_deletion,chg_codon,pos_sequence,type_nt,info)
                    
                    print("a_bis",ref_upper,alt_upper)
                    
                    
                    print("b",sequence[first_pos-1])
                    
                    print("c",shift)
                    
                    print("d",premier,dernier)
                    
                    print("e",suivant)
            
                    
                    
                print("f",p,d,s)
                
                
                print("\n")
    
    
    if chg_codon!="":
    
    
        if "-"  in str(pos_sequence) and "-" not in nt and only_gap(separator(chg_codon)[0]):
            a=0
            
            ##########################
            #preparation of variables#
            ##########################
            
            ref,alt=separator(chg_codon)
            
            #compute the lg of event according to the position SNP reference
            lg_even=lg_SNP_compute(pos_sequence)
            
            #first position in the reference codon  123-125 <- 123
            first_pos=SNPposition_to_consider(pos_sequence)
            
            first_decoup_seq=first_pos+lg_even-1
            
            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
            
            nt_alt=pos_ajust_alt(alt,nt)
    
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################  
            
            seq_alt=sequence[first_pos:first_pos+len(alt)]
            
    
            if alt==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0
    
    
                
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
                
                ev+="&nothing_insert"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
                
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1:
                #ok definitif
    
                ev+="&deletion_ajustee"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos
                lg_deletion=len(seq_alt)
    
                
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0:
                #ok definitif
                ev+="&new_insert_spe"
                even_choix="insert"           
                pos_sequence=first_pos
    
                #warning fir insertion, it'possible to find nucleotides in too much in nt compared to alt
                #we need remove nucleotides in too much
                
                if alt in nt:
                    
                    nt=alt
                    
                else:
    
                    even_choix="upgrade"
                    ev="upgrade"
             
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
                            
                ev+="&nothing_insert_empty"
                even_choix="nothing_insert"
                pos_sequence=first_pos
                
        elif "-"  in str(pos_sequence) and "-" not in nt and not only_gap(separator(chg_codon)[0]):
            a=0
            
            
            ##########################
            #preparation of variables#
            ##########################
            
            ref,alt=separator(chg_codon)
            
    
            #compute the lg of event in terms of position of SNP
            lg_even=lg_SNP_compute(pos_sequence)
            
    
            #first position in the reference codon  123-125 <- 123
            first_pos=SNPposition_to_consider(pos_sequence)
            
            first_decoup_seq=first_pos+lg_even-1
            
            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
    
            pos_ajust,nt_alt=pos_ajust_alt_2(alt,nt)
            
    
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################  
            
            seq_alt=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(alt)-1-pos_ajust]
            
    
            if alt.upper()==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0   
    
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
    
                ev+="&nothing_insert_2"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
             
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1 :
                
                ev+="&deletion_ajustee_2"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos+pos_ajust-1
    
                lg_deletion=count_maj(alt)
                
    
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0 :
               
                ev+="&new_insert_spe_2"
                even_choix="insert"
    
                pos_sequence=first_pos
    
                if enlever_minuscule(alt) in nt:
                    
                    nt=enlever_minuscule(alt)
                    
                else:
                    
                    even_choix="upgrade"
                    ev="upgrade"
            
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
    
    
                ev+="&nothing_insert_empty_2"
                even_choix="nothing_insert"
                pos_sequence=first_pos
                
                
        elif "-" not in str(pos_sequence) and "-" not in nt and "-" not in chg_codon:
            a=0
            
            ###########################
            #preparation des variables#
            ###########################
            
            ref,alt=separator(chg_codon)
            
            #compute the lg of event in terms of position of SNP
            lg_even=len(ref)-len(alt)
            
    
            #first position in the reference codon  123-125 <- 123
            first_pos=pos_sequence
            
            first_decoup_seq=first_pos+lg_even-1
            
            ########################
            #compare the nt and alt#
            #using alt because it is a insertion#
            #####################################
    
            pos_ajust,nt_alt=pos_ajust_alt_2(alt,nt)
            
    
            #################################################################
            #comparison ot the nt of the reference sequence with ref and alt#
            #################################################################    
            
            seq_alt=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(alt)-1-pos_ajust]
            
    
            if alt.upper()==seq_alt:
                
                choice_alt_seq=1
                
            else:
                
                choice_alt_seq=0   
    
            #################
            #Choose of event#
            #################
            
            #the nt is a insertion and it's exist in the sequence, nothing to do
            if nt_alt==1 and choice_alt_seq==1:
    
                ev+="&nothing_insert_3"
                even_choix="nothing_insert"
                pos_sequence=first_pos
    
             
            #the nt is not a insestion but it's exists in the sequence, it must be removed
            elif nt_alt==0 and choice_alt_seq==1 :
                
                ev+="&deletion_ajustee_3"
                even_choix="deletion_ajustee"
                pos_sequence=first_pos+pos_ajust-1
    
                lg_deletion=count_maj(alt)
                
    
            #the nt is a insertion but the insertion doesn't exits in the sequence, we must insert them.
            elif nt_alt==1 and choice_alt_seq==0 :
               
                ev+="&new_insert_spe_3"
                even_choix="insert"
    
                pos_sequence=first_pos
    
                if enlever_minuscule(alt) in nt:
                    
                    nt=enlever_minuscule(alt)
                    
                else:
                    
                    even_choix="upgrade"
                    ev="upgrade"
            
            #the read is not a insertion and it's no existing in the sequence
            elif nt_alt==0 and choice_alt_seq==0:
    
    
                ev+="&nothing_insert_empty_3"
                even_choix="nothing_insert"
                pos_sequence=first_pos

        else : 
            
            even_choix="upgrade"
            ev="upgrade"
            
                
    elif chg_codon=="" and len(nt)==1:
        a=0
        
        if "-" in str(pos_sequence):
        
            lg_even=lg_SNP_compute(pos_sequence)
            first_pos=SNPposition_to_consider(pos_sequence)
            
        else:
            
            lg_even=1
            first_pos=int(pos_sequence)

        nt_suivant_even=sequence[first_pos-1:first_pos-1+lg_even]
        
        nt_courant=sequence[first_pos-1]
        
        nt_suivant=sequence[first_pos]
        
        
        #si nt==nt_courant, il n'y a rien a faire car la variant est deja considere
        
        if nt==nt_courant:
            
            ev+="&chg_codon_empty_nothing"
            even_choix="nothing_delet"
            pos_sequence=first_pos
            
        #si nt==nt_suivant, le nt est le suivant donc il faut faire une deletion
            
        elif nt==nt_suivant:
            
            ev+="&deletion_chg_codon_empty"
            even_choix="deletion_ajustee"
            pos_sequence=first_pos
            lg_deletion=len(nt)
            
        else : 
    
            even_choix="upgrade"
            ev="upgrade"   
            
        
            
    else : 

        even_choix="upgrade"
        ev="upgrade"   
        
    return even_choix,pos_sequence,nt,lg_deletion,ev


def choix_particulier_deletion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info):
    ''' documentation in preparation '''  

    if "-"  in str(pos_sequence) and "-" not in nt and only_gap(separator(chg_codon)[1]):
        

        ##########################
        #preparation of variables#
        ##########################
        
        ref,alt=separator(chg_codon)
        
        #compute the lg of event in terms of position of SNP
        lg_even=lg_SNP_compute(pos_sequence)
        
        #premiere position dans le codon de reference  123-125 <- 123
        first_pos=SNPposition_to_consider(pos_sequence)
        
        first_decoup_seq=first_pos+lg_even-1
        
        ########################
        #compare the nt and alt#
        #using alt because it is a insertion#
        #####################################
        
        nt_ref=pos_ajust_alt(ref,nt)
     
        #################################################################
        #comparison ot the nt of the reference sequence with ref and alt#
        #################################################################  

        seq_ref=sequence[first_pos-1:first_pos+len(ref)-1]


        if ref==seq_ref:
            
            choice_ref_seq=1
            
        else:
            
            choice_ref_seq=0
            
        #################
        #Choose of event#
        #################
        

        
        #the nt is a deletion and it's exist in the sequence, nothing to do
        if nt_ref==1 and choice_ref_seq==1:
            

            ev+="&nothing_delet"
            even_choix="nothing_delet"
            pos_sequence=first_pos
            
        #the nt is not a deletion but it's exists in the sequence, it must be removed
        elif nt_ref==0 and choice_ref_seq==1:

            ev+="&deletion_spe"
            even_choix="deletion_ajustee"
            pos_sequence=first_pos
            lg_deletion=len(seq_ref)
            
            
        else:
            
            

            even_choix="upgrade_deletion_1"
            ev="upgrade"
       
    elif "-"  in str(pos_sequence) and "-" not in nt and not only_gap(separator(chg_codon)[1]):  
        
         
        ##########################
        #preparation of variables#
        ##########################
        
        ref,alt=separator(chg_codon)
        
        #compute the lg of event in terms of position of SNP
        lg_even=lg_SNP_compute(pos_sequence)
        

        #first position in the reference codon  123-125 <- 123
        first_pos=SNPposition_to_consider(pos_sequence)
        
        first_decoup_seq=first_pos+lg_even-1
        
        ########################
        #compare the nt and alt#
        #using alt because it is a insertion#
        #####################################

        pos_ajust,nt_ref=pos_ajust_alt_2(ref,nt)
        

        #################################################################
        #comparison ot the nt of the reference sequence with ref and alt#
        #################################################################   
        
        seq_ref=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(ref)-1-pos_ajust]
        
        
        if ref.upper()==seq_ref:
            
            choice_ref_seq=1
            
        else:
            
            choice_ref_seq=0 
            
            
        #the nt is a deletion and it's exist in the sequence, nothing to do
        if nt_ref==1 and choice_ref_seq==1:

            ev+="&nothing_delet_2"
            even_choix="nothing_delet"
            pos_sequence=first_pos
            
        #the nt is not a deletion but it's exists in the sequence, it must be removed
        elif nt_ref==0 and choice_ref_seq==1 :
            
            ev+="&deletion_spe_2"
            even_choix="deletion_ajustee"
            pos_sequence=first_pos
            lg_deletion=count_maj(ref)
            
            
        #the nt is a deletion but the insertion doesn't exits in the sequence, we must insert them
        elif nt_ref==1 and choice_ref_seq==0 :

            ev+="&insert_de_deletion_2"
            even_choix="insert"
            pos_sequence=first_pos+1
            nt=enlever_minuscule(ref)
            
        #the read is not a deletion and it's no existing in the sequence
        elif nt_ref==0 and choice_ref_seq==0 :   
    
            ev+="&nothing_delet_empty_2"
            even_choix="nothing_delet"
            pos_sequence=first_pos
            
            
    elif "-" not in str(pos_sequence) and "-" not in nt and "-" not in chg_codon:

        ##########################
        #preparation of variables#
        ##########################
        
        ref,alt=separator(chg_codon)
        
        lg_even=len(ref)-len(alt)
    
        #first position in the reference codon  123-125 <- 123
        first_pos=pos_sequence
        
        first_decoup_seq=first_pos+lg_even-1
        
        ########################
        #compare the nt and alt#
        #using alt because it is a insertion#
        #####################################

        pos_ajust,nt_ref=pos_ajust_alt_2(ref,nt)
        
        #################################################################
        #comparison ot the nt of the reference sequence with ref and alt#
        #################################################################   
        
        seq_ref=sequence[first_decoup_seq-1-pos_ajust:first_decoup_seq+len(ref)-1-pos_ajust]
               
        if ref.upper()==seq_ref:
            
            choice_ref_seq=1
            
        else:
            
            choice_ref_seq=0 
                       
        #the nt is a deletion and it's exist in the sequence, nothing to do
        if nt_ref==1 and choice_ref_seq==1:

            ev+="&nothing_delet_3"
            even_choix="nothing_delet"
            pos_sequence=first_pos
            
        #the nt is not a deletion but it's exists in the sequence, it must be removed
        elif nt_ref==0 and choice_ref_seq==1 :
            
            ev+="&deletion_spe_3"
            even_choix="deletion_ajustee"
            pos_sequence=first_pos
            lg_deletion=count_maj(ref)
                       
        #the nt is a deletion but the insertion doesn't exits in the sequence, we must insert them
        elif nt_ref==1 and choice_ref_seq==0 :
        
            ev+="&insert_de_deletion_3"
            even_choix="insert"
            pos_sequence=first_pos+1
            nt=enlever_minuscule(ref)
            
        #the read is not a deletion and it's no existing in the sequence
        elif nt_ref==0 and choice_ref_seq==0 :   
    
            ev+="&nothing_delet_empty_3"
            even_choix="nothing_delet"
            pos_sequence=first_pos

    else:
         
        even_choix="upgrade_deletion_2"
        ev="upgrade"
        
    return even_choix,pos_sequence,nt,lg_deletion,ev


def choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info):
    """
    choix_particulier_frameshift
    """
    ''' documentation in preparation '''    
    if chg_codon!="":
    
        ref,alt=separator(chg_codon)
        
        count_ref=count_maj(ref)
        
        count_alt=count_maj(alt)
    
        
        if "-" in nt:
        
            if count_ref > count_alt:
                
                even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_deletion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
                
            elif count_ref < count_alt:
                                
                even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
                
        elif "-" not in nt:
            
    
            
            if count_ref > count_alt:
                even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_deletion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
    
            
            elif count_ref < count_alt:
                                
                even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
    
        else:
    
        
            even_choix="upgrade_frameshift"
            ev="upgrade"
            
    else:
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
    
    return even_choix,pos_sequence,nt,lg_deletion,ev

#################################
#functions of application events#
#################################
    
def insert_seq(sequence,nt,pos_sequence,ev,info):
    #return a sequence with a insertion
    #insert the element just after the position given
    #if "-" in position, give only first position
    ''' documentation in preparation '''
    return sequence[:pos_sequence]+nt+sequence[pos_sequence:]

def as_stop(seq):
    #return 1 if the given sequence is a codon stop
    ''' documentation in preparation '''   
    if seq=="TAA" or seq=="TGA" or seq=="TAG":
        
        retour=1
        
    else:
        
        retour=0
    
    return retour

def remplacement(sequence,D,pos_SNP):
    """
    replace a nt in the given sequence by another nt given
    """
    ''' documentation in preparation ''' 
    return sequence[:pos_SNP-1]+D+sequence[pos_SNP+len(D)-1:]

def deletion(sequence,pos_sequence,lg_deletion):
    #remove n nucleotide at the given position
    #if "-" in position, give only first position
    ''' documentation in preparation ''' 
    return sequence[:pos_sequence]+sequence[pos_sequence+lg_deletion:]

#######################
#Specific applications#
#######################
    
def find_first_start(sequence):
    """
    return the sequence in input with all nt before the first start codon cutted
    """
    ''' documentation in preparation '''
    start=Seq(sequence).find("ATG")

    return sequence[start:],start

########################
#functions cut sequence#
########################
    
def decoupage_stop(pos_snp,sequence,lgsequence,profondeur,info):
    ''' documentation in preparation '''
    #Can cut the sequence in downstream of stop lost until the next stop

    #Compute upstream and downstream position
    amont=pos_snp-profondeur-position_SNP(pos_snp)
    aval=len(sequence)
    

    return sequence[amont-1:aval],(pos_snp + ( position_SNP_inverse(pos_snp) % 3 ))/3
    
def position_SNP(pos_SNP):
    #return a position of SNP in his codon
    ''' documentation in preparation '''    

    if (pos_SNP) % 3 == 1:
        
        return 0
    
    elif (pos_SNP) % 3 == 2:
        
        return 1
    
    elif (pos_SNP) % 3 == 0:
        
        return 2
    
def position_SNP_inverse(pos_SNP):
    #return a position of SNP in his codon
    ''' documentation in preparation '''
    if (pos_SNP) % 3 == 1:
        
        return 2
    
    elif (pos_SNP) % 3 == 2:
        
        return 1
    
    elif (pos_SNP) % 3 == 0:
        
        return 0
    
def decoupage_insert_v2(pos_SNP,sequence,lg,profondeur,info,ev,nt):
    ''' documentation in preparation '''
    #cut a sequence according to both position, begin and end of insert 

    ref_c=enlever_minuscule(info.split("|")[8].split("/")[0])
    alt_c=enlever_minuscule(info.split("|")[8].split("/")[1])
    
    lg_nt=max(len(ref_c),len(alt_c))

    pos_peptide=0

    #Compute positions of begin and end of insertion

    pos_begin=pos_SNP+1
    pos_end=pos_SNP+lg_nt
    
    
    #Compute codon positions of begin and end of insertion
    
    pos_codon_begin=position_SNP(pos_begin)
    pos_codon_end=position_SNP_inverse(pos_end)
    
    
    if pos_begin<=profondeur:
        
        ajust_amont=1
        
    else:
        
        ajust_amont=pos_begin-profondeur-pos_codon_begin
        
     
    if pos_end>=len(sequence):
        
        ajust_aval=len(sequence)
        
    else:
        
        ajust_aval=pos_end+pos_codon_end+profondeur

    sequence=sequence[ajust_amont-1:ajust_aval]
    pos_peptide=(pos_SNP + ( position_SNP_inverse(pos_SNP) % 3 ))/3
        

    return sequence,pos_peptide

def decoupage_v2(pos_SNP,sequence,lg,profondeur,info):
    #cut a sequence around a given position 
    ''' documentation in preparation '''    
    pos_protein=0
    #position du nt dans le codon
    pos_codon=position_SNP(pos_SNP)
    

    if pos_SNP % 3 == 1:

        if pos_SNP<=profondeur:
            
            ajust_amont=1

        else:
        
            ajust_amont=pos_SNP-profondeur-pos_codon
        
        if pos_SNP+position_SNP_inverse(pos_SNP)+profondeur>=len(sequence):
                     
            ajust_aval=len(sequence)
            
        else:
        
            ajust_aval=pos_SNP+position_SNP_inverse(pos_SNP)+profondeur

        sequence=sequence[ajust_amont-1:ajust_aval]


       
    elif pos_SNP % 3 == 2:
        
        if pos_SNP<=profondeur:
            
            ajust_amont=1

        else:
        
            ajust_amont=pos_SNP-pos_codon-profondeur

        if pos_SNP+position_SNP_inverse(pos_SNP)+profondeur>=len(sequence):
                     
            ajust_aval=len(sequence)
            
        else:
        
            ajust_aval=pos_SNP+position_SNP_inverse(pos_SNP)+profondeur
    
        sequence=sequence[ajust_amont-1:ajust_aval]



    elif pos_SNP % 3 == 0:

        if pos_SNP<=profondeur:
            
            ajust_amont=1

        else:
        
            ajust_amont=pos_SNP-pos_codon-profondeur
                     
        if pos_SNP+position_SNP_inverse(pos_SNP)+profondeur>=len(sequence):
                     
            ajust_aval=len(sequence)
            
        else:
        
            ajust_aval=pos_SNP+profondeur+position_SNP_inverse(pos_SNP)
            

        sequence=sequence[ajust_amont-1:ajust_aval]

        
    if ajust_amont==1:
        
        pos_protein=(pos_SNP+position_SNP_inverse(pos_SNP))/3
        
        
    else:    
    
        pos_protein=(pos_SNP+position_SNP_inverse(pos_SNP)-ajust_amont-position_SNP_inverse(ajust_amont)+3)/3

    
    return sequence,pos_protein


################
#Main functions#
################
    
def choix_particulier_v2(seq_stop_upgrade,nt,lg,sequence,consequence,pos_sequence,chg_codon,gene,transcrit,ev,info,seq_stop):
    """
    choix_particulier_v2
    """
    ''' documentation in preparation '''     

    even_choix="None"
    lg_deletion=0
    add_seq=""

    
    if "missense_variant" in consequence and len(consequence)==1:
        
        ev="missense_spe"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_missense(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence)
        
    elif "inframe_insertion" in consequence and len(consequence)==1:
   
        ev="insertion"
                
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
    elif "inframe_insertion" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
        
        ev="insertion"
                
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
    elif "inframe_deletion" in consequence and len(consequence)==1:
        
        ev="deletion"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_deletion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
    
    elif "frameshift_variant" in consequence and len(consequence)==1:

        ev="deletion"

        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
    elif "frameshift_variant" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
        
        ev="frameshift"

        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "frameshift_variant" in consequence and len(consequence)==1 and "-" not in str(pos_sequence):
        
        ev="frameshift"

        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
    elif "stop_gained" in consequence and len(consequence)==1:
        
        ev="stop_gained_spe"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "stop_gained" in consequence and "frameshift_variant" in consequence and len(consequence)==2:
        
        ev="frameshift_stop_gained"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "stop_gained" in consequence and "inframe_insertion" in consequence and len(consequence)==2:
        
        ev="stop_gained_inframe_insertion"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)

    elif "stop_lost" in consequence: 
        
            ev="frameshift_variant_stop_lost"
        
        
            if "-" in str(pos_sequence):
    
                pos_sequence=SNPposition_to_consider(pos_sequence)
                
            sequence=remplacement(sequence,nt,pos_sequence)
    
            if as_stop(sequence[-3:]):
                    
                even_choix="Nothing"

                    
            else:
                    
                #find the next stop
                pos_next_stop,seq_stop_upgrade=det_stop(seq_stop)
                #cut the sequence stop according to the next stop
                add_seq=cut_seq(seq_stop,pos_next_stop)
        
                even_choix="ajout_sequence"


    elif "protein_altering_variant" in consequence : 
        
        ev="protein_altering_variant"
    
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "inframe_deletion" in consequence and "start_lost" and len(consequence)==2:
        
        ev="start_lost_and_inframe_deletion"
    
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_deletion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "frameshift_variant" in consequence and "stop_gained" in consequence and len(consequence)==3:
        
        ev="frameshift_variant_and_stop_retained_variant"
    
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
    elif "start_retained_variant" in consequence:   

        ev="start_retained_variant"
        
        even_choix="start_retained_variant"
                
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_insertion(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
                
            
    elif "start_lost" in consequence and "frameshift_variant" in consequence and len(consequence)==2:
        
        ev="frameshift_variant&lose_start"
        
        even_choix="lose_start_particulier"
    
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "frameshift_variant" in consequence and "feature_truncation" in consequence and len(consequence)==2 and chg_codon=="":
                
        ev="frameshift_variant&feature_truncation_without_chg_codon"
        
        even_choix="frameshift_variant"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
        
    elif "frameshift_variant" in consequence and "feature_elongation" in consequence and len(consequence)==2 and chg_codon=="":
                
        ev="frameshift_variant&feature_elongation_without_chg_codon"
        
        even_choix="frameshift_variant"
        
        even_choix,pos_sequence,nt,lg_deletion,ev=choix_particulier_frameshift(even_choix,pos_sequence,nt,lg_deletion,ev,chg_codon,sequence,info)
        
 
    else:
        
        even_choix="upgrade_choix_particulier"
        ev="upgrade"
        
    
    
    return even_choix,pos_sequence,nt,lg_deletion,ev,add_seq,seq_stop_upgrade



def construction_peptide_v2(nt,lg,sequence,consequence,pos_sequence,chg_codon,seq_stop,gene,transcrit,info):
    """
    construction_peptide_v2
    """
    ''' documentation in preparation '''    

    #function none efficient for  "splice_acceptor_variant" and "splice_donor_variant"

    #even variable
    even_choix=""
    
    #SNP position into peptide
    pos_SNP_seq=0
    
    #lg_deletion 
    lg_deletion=0
    
    #seq stop
    add_seq=""
    
    seq_stop_upgrade=0

    ######################
    #choix de l'evenement#
    ######################

    #output variable of even
    ev="chgt_simple"
    

    
    if len(chg_codon)==7 and "-"  not in nt and len(nt)==1 and "-" not in str(pos_sequence):
        

        #only missense_variant
        if "missense_variant" in consequence and len(consequence)==1: 
            
            even_choix="chgt_simple"
                    
                
        elif "missense_variant" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
                
            even_choix="chgt_simple"
                
                
        elif "initiator_codon_variant" in consequence and len(consequence)==1:
                
            even_choix="chgt_simple"
                
        elif ("stop_lost" in consequence and len(nt)==1) or ("frameshift_variant" in consequence and "stop_lost" in consequence and len(consequence)==2):    
                
            sequence=remplacement(sequence,nt,pos_sequence)
    
            if as_stop(sequence[-3:]):
                    
                even_choix="Nothing"

                    
            else:
                    
                #find the next stop
                pos_next_stop,seq_stop_upgrade=det_stop(seq_stop)
                
                
                #cut the stp sequence at next stop
                add_seq=cut_seq(seq_stop,pos_next_stop)

        
                even_choix="ajout_sequence"    
                
        elif 'initiator_codon_variant' in consequence and 'splice_region_variant' in consequence and len(consequence)==2:
                
            even_choix="chgt_simple"

            
        elif "initiator_codon_variant" in consequence and "NMD_transcript_variant" in consequence and len(consequence)==2:
            
            even_choix="chgt_simple"
            
        elif "stop_gained" in consequence and len(consequence)==1:
                
            even_choix="chgt_simple_gain_stop"
                
        elif "stop_gained" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
                
            even_choix="chgt_simple_gain_stop"
                   
        elif "stop_gained" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
                
            even_choix="chgt_simple_gain_stop"
            
        elif "start_lost" in consequence and len(consequence)==1:
            
            
            even_choix="chg_simple_lose_start"
            ev="start_lose"
            
        elif "start_lost" in consequence and "splice_region_variant" in consequence and len(consequence)==2:
            
            even_choix="chg_simple_lose_start"
            ev="start_lose"
            

        else:
            
            #select NMD variant
            if "NMD_transcript_variant" in consequence and "missense_variant" in consequence and len(consequence)==2:
                
                even_choix="chgt_simple"
                
            elif "stop_gained" in consequence and "NMD_transcript_variant" in consequence and len(consequence)==2:
                
                even_choix="chgt_simple_gain_stop"
                
            else:
                
                even_choix="upgrade_choix_classique"
                ev="upgrade"
    
                pos_SNP_seq=0

    #Particular case treat in another function
    else:
       
        ev="cas_specifique"
        
        even_choix,pos_sequence,nt,lg_deletion,ev,add_seq,seq_stop_upgrade=choix_particulier_v2(seq_stop_upgrade,nt,lg,sequence,consequence,pos_sequence,chg_codon,gene,transcrit,ev,info,seq_stop)
        
            
    #############
    #apply event#
    #############
    
    if even_choix=="chgt_simple":
        
        #change the nt
        sequence=remplacement(sequence,nt,pos_sequence)
        
    elif even_choix=="chg_simple_lose_start":
        
        #change the nt
        sequence=remplacement(sequence,nt,pos_sequence)
              
    elif even_choix=="ajout_sequence":    

        #add sequecnce after stop lost
        sequence=remplacement(sequence+add_seq,nt,pos_sequence) 


    elif even_choix=="chgt_simple_gain_stop":
            
        sequence=remplacement(sequence,nt,pos_sequence)
                
    elif even_choix=="insert":
        
        sequence=insert_seq(sequence,nt,pos_sequence,ev,info)
   
    elif even_choix=="deletion":

        sequence=deletion(sequence,pos_sequence,lg_deletion)
  
    elif even_choix=="deletion_ajustee":

        sequence=deletion(sequence,pos_sequence,lg_deletion)
        

    elif even_choix=="None":
            
        sequence=""
        pos_SNP_seq=0

    #######################
    #Specific applications#
    #######################
    
    if "start_retained_variant" in ev:
        
        sequence=find_first_start(sequence)[0]
        

    ######################
    #Cutting the sequence# 
    ######################

    #decision_aval=d_aval(sequence,first_sequence,pos_sequence)
    
        
    if even_choix=="chgt_simple":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="chg_simple_lose_start":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="ajout_sequence":
        
        sequence,pos_SNP_seq=decoupage_stop(pos_sequence,sequence,len(sequence),(lg-1)*3,info)

    elif even_choix=="frameshift_variant":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)   
    
    elif even_choix=="nothing_insert":
        
        sequence,pos_SNP_seq=decoupage_insert_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info,ev,nt)

    elif even_choix=="chgt_simple_gain_stop":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)

    elif even_choix=="insert":
        
        sequence,pos_SNP_seq=decoupage_insert_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info,ev,nt)

    elif even_choix=="deletion_insertion":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="deletion":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="Nothing":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="deletion_ajustee":

        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="nothing_delet":
        
        sequence,pos_SNP_seq=decoupage_v2(pos_sequence,sequence,len(sequence),(lg-1)*3,info)
        
    elif even_choix=="None":
            
        sequence=""
        pos_SNP_seq=0
        
        
    ###########
    #translate#
    ###########     
    
        
    if "lose_start" in even_choix:
    
        pos_start=find_first_start(sequence)[1]
        
        if pos_start==0:
            
            if "N" not in sequence:
            
                sequence=Seq(sequence, IUPAC.unambiguous_dna).translate(to_stop=True)      
                        
            elif "N" in sequence:
                        
                sequence=Seq(sequence, IUPAC.ambiguous_dna).translate(to_stop=True)     
        
            elif even_choix=="None":
                  
                sequence=""
                pos_SNP_seq=0
                
        else:
            
            sequence=""
            pos_SNP_seq=0
            

    else:
           
    
        if "N" not in sequence:

                
            sequence=Seq(sequence, IUPAC.unambiguous_dna).translate(to_stop=True)    
            

        elif "N" in sequence:

            sequence=Seq(sequence, IUPAC.ambiguous_dna).translate(to_stop=True)  

        elif even_choix=="None":
              
            sequence=""
            pos_SNP_seq=0

    

    return sequence,int(pos_SNP_seq),ev,seq_stop_upgrade,even_choix

