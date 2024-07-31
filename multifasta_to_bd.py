#!/usr/bin/env python3

# Ce programme a été créé par Yascim Kamel dans le cadre de son stage de Master 2 au cours du projet ANR DICWOC porté
# par Aurélie Célérier et Sylvia Campagna, et modifié par Ambre Petit dans le cadre de son stage de Master 1 à l'Université de Montpellier. Il permet de lancer les programmes d'implémentation de la base de données associée au projet.

import Bio.pairwise2 as pairwise2
import re
import sys
import psycopg2
import Bio.SeqIO as SeqIO


########### Fonctions d'appel de toutes les autres fonctions  ###########

def multi_fasta_to_bd(multi_fasta, connection, experience_ID, espece, genre, BD, assemblie_ID):
    conversion_multifasta(multi_fasta)
    insertion_organism(connection, genre, espece)
    insertion_assemblie(connection, genre, espece, BD, assemblie_ID)
    insertion_gene(tableau_gene, connection, experience_ID, espece,assemblie_ID)


########### Conversion du fichier multifasta en un tableau à deux dimensions  ###########

def conversion_multifasta(multi_fasta):
    global tableau_gene
    tableau_gene=[]

    for record in SeqIO.parse(multi_fasta,"fasta"): #Utilisation de SeqIO de Biopython pour parcourir les séquences fasta contenues dans le multifasta
        header=str(record.id) #Récupération du header de la séquence fasta en chaîne de caractères.
        split_header=header.split('|') #Récupération des informations contenues dans le header
        nom=split_header[0]
        famille=split_header[1]
        if len(split_header)==3 : #L'état de pseudogène apparait dans le header contrairement à l'état fonctionnel
            etat="pseudogène"
        else: #Si pas d'information dans de header alors le gène est par défaut fonctionnel
            etat="fonctionnel"
        position= nom.split(':')
        liste_position= position[1].split("-")
        start=liste_position[0] #Position de début du gène OR
        end=liste_position[1] #Position de fin du gène OR

        gene=["",nom,'gène OR',famille,etat,int(start),int(end),str(record.seq),"" ]
        tableau_gene.append(gene) #Stockage des différents gènes dans un tableau à deux dimensions
        
        
########### Remplissage de la table "organism" ###########

def insertion_organism(connection, genre, espece):
    cur=connection.cursor()
    query_organism="INSERT INTO public.organism(\"Genre\",\"Espèce\") SELECT \'"+str(genre)+"\',\'"+str(espece)+"\' WHERE NOT EXISTS (SELECT 1 FROM public.organism WHERE \"Genre\" = \'"+str(genre)+"\' AND \"Espèce\" = \'"+str(espece)+"\');" #Le remplissage se fait seulement si les informations sur l'organisme n'est pas déjà présent dans la table
    cur.execute(query_organism)
    connection.commit()
    # Affichage du message si l'insertion a réussi
    if cur.rowcount > 0:
        print("\n Insertion des données de", genre, espece, "dans la table organism ...")
    
    
########### Remplissage de la table "assemblie"  ###########

def insertion_assemblie(connection, genre, espece, BD, ID_A):
    cur=connection.cursor()
    query_assemblie="INSERT INTO public.assemblie(\"ID\", \"Genre\",\"Espèce\",\"Base de donnee\") SELECT \'"+str(ID_A)+"\',\'"+str(genre)+"\',\'"+str(espece)+"\',\'"+str(BD)+"\' WHERE NOT EXISTS (SELECT 1 FROM public.assemblie WHERE \"ID\" = '"+str(ID_A)+"');" #Le remplissage se fait seulement si les informations sur l'assemblage n'est pas déjà présent dans la table
    cur.execute(query_assemblie)
    connection.commit()
    # Affichage du message si l'insertion a réussi
    if cur.rowcount > 0:
        print("\n Insertion des données de", genre, espece, "dans la table assemblie ...")
    

########### Remplissage de la table "gene" en fonction des paramètres d'alignement de séquences (calcul de similarité)  ###########

def insertion_gene(tmf,connection,ID_E,esp, ID_A):
    print ("\n Debut de l'insertion des données dans les tables link et gene ...")
    #### Définition des valeurs utilisées dans l'alignement
    match=2
    missmatch=-2
    gap_int=-1
    gap_ext=0

    cur=connection.cursor()
    
    gene_ref="SELECT distinct link.\"ID_gene\", gene.\"Nom\", gene.\"Séquence\" FROM link,gene WHERE link.\"ID_gene\"= gene.\"ID\" AND link.\"ID_assemblie\" IN (SELECT assemblie.\"ID\" FROM assemblie WHERE \"Espèce\"='"+espece+"') AND link.\"ID_gene\" in (SELECT gene.\"ID\" FROM gene WHERE gene.\"ID\" = gene.\"Référence\");" #Requete qui permet de récupérer seulement les gènes référents d'un meme assemblage pour une espèce donnée

    cur.execute(gene_ref)
    tc=cur.fetchall()
    IDg=0
    cur.execute("SELECT max(\"ID\") FROM gene")
    ID_max=cur.fetchall() #Récupération du dernier ID de gène généré
    
    # Vérifier si la base de données est vide
    cur.execute("SELECT COUNT(*) FROM gene")
    count = cur.fetchone()[0]
    
    if count == 0:
    #Si base de données vide, insertion des gènes sans comparaison puisqu'il n'y a pas de gènes référents
        comparaison_insertion(tmf, connection, cur, IDg, tc, esp, ID_A, ID_E, ID_max)
    else :
    #Insertion des gènes avec étape préliminaire de comparaison avec les gènes référents
        IDg=ID_max[0][0]
        comparaison_insertion(tmf, connection, cur, IDg, tc, esp, ID_A, ID_E, ID_max)

    
def comparaison_insertion(tmf,connection, cur, IDg, tc, esp, ID_A, ID_E, ID_max): #Permet de comparer les séquences du gène d'interet et des références afin d'adapter l'insertion à la base de données
    for gt in tmf: # Parcours de la liste des gènes trouvés par le pipeline
        print("\nTraitement du gène : ", gt[1], " Assemblage : "+ ID_A + "\n")
        print('\n_________________________________________________________________________\n')
        gr_id_max=""
        gr_score_max=0
        sim_max=0
        nb_match=0
        IDg=IDg+1
        list_query_gene="INSERT INTO public.gene VALUES "
        list_query_link="INSERT INTO public.link VALUES "
        
        if tc: # Comparaison avec les gènes de références seulement s'ils existent
            for gr in tc: # Parcours des gènes de références
                lgr=len(gr[2])
                lgt=len(gt[-2])
                print("Comparaison entre :", gt[1] , " et ", gr[1] , "…")
                score_match=align(gt[-2], gr[2] , -1 ,-1) # Appel de la fonction align pour effectuer l'alignement de la sequence d'interet et des sequences de références
                print("Score d alignement : ", score_match[0], "Nombre de match : ", score_match[1])
                if score_match[0] > gr_score_max: #Boucle qui permet de calculer la similarité seulement si le score d'alignement est plus grand que le dernier score d'alignement maximum
                    gr_id_max=gr[0]
                    gr_score_max=score_match[0]
                    nb_match=score_match[1]
                    if lgt<lgr:
                        sim_max=nb_match/lgt
                        print("Similarité entre les séquences : ", sim_max)
                    else :
                        sim_max=nb_match/lgr
                        print("Similarité entre les séquences : ", sim_max)
        
        
        if sim_max>0.90: # Si le gène d'interet est similaire à 90 % à un gène de référence, alors ce gène devient la référence du gène d'intéret et la séquence n'est pas rentrée dans la base de donnée
                list_query_gene+="("+str(IDg)+",'"+ str(gt[1])+"','gene_OR','"+ str(gt[3])+"','"+ str(gt[4])+"'," +str(gt[5])+","+ str(gt[6])+",'Se référer au gène de référence',"+str(gr_id_max)+");"
        else : # Par défaut sim_max vaut 0 donc si les gènes de références n'existent pas, les gènes sont rentrés sans faire référence aux gènes de références
                list_query_gene+="("+str(IDg)+",'"+ str(gt[1])+"','gene_OR','"+ str(gt[3])+"','"+ str(gt[4])+"'," +str(gt[5])+","+ str(gt[6])+",'"+str(gt[7])+"',"+str(IDg)+");"
        
        list_query_link+="('"+str(ID_A)+"',"+str(ID_E)+","+str(IDg)+");" #Remplissage de la table Link pour chaque gène de l'assemblage
        
        cur.execute(list_query_gene)
        connection.commit()

        cur.execute(list_query_link)
        connection.commit()

        print(gt[1], " Traité")
        print('\n_________________________________________________________________________\n')

    print("\nAssemblie " + assemblie_ID +" entièrement traité.")


########### Fonctions pour l'alignement des deux séquences  ###########

def align(seq1, seq2, open_gap, extend_gap):
    score_alignment = pairwise2.align.globalxs(seq1, seq2,open_gap,extend_gap,penalize_end_gaps=(False,False))
    
    return number_match(score_alignment[0][0], score_alignment[0][1], score_alignment[0][2])

def number_match(align1, align2, score_alignment): #Compte le nombre de match entre deux séquences
    s = [score_alignment]
    match=0
    i=0
    if len(align1)==len(align2):
        while i<len(align1):
            if align1[i]==align2[i] and (align1[i]!="-" or align2[i]!="-"):
                match+=1
            i+=1
        s.append(match)
    else:
        print("Erreur : les alignements passés en argument ne sont pas de la même longueur ")
    return s



###########  Main  ###########

connection = psycopg2.connect("dbname='CeGeC' user=Ambre host='localhost' port='5432'")

## Récupération des paramètres lors de l'appel du script :
nom_complet=sys.argv[1] #Récupération du nom de l'organisme
BD=sys.argv[2] #Récupération du nom de la base de données
esp_genre=nom_complet.split("_") #Découpe du nom de l'assemblage en deux parties
espece=esp_genre[1]
genre=esp_genre[0]
assemblie_ID= nom_complet +'_'+ BD #Creation de l'ID de l'assemblage en fonction du nom de l'organisme et de la base de données
multi_fasta=sys.argv[3] #Fichier fasta des gènes OR de l'organisme d'intéret
experience_ID=sys.argv[4] #ID de l'experience généré avec le Snakemake

multi_fasta_to_bd(multi_fasta,connection,experience_ID,espece, genre, BD, assemblie_ID)

connection.close()
