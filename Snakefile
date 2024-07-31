########### Snakefile créé par Yascim Kamel dans le cadre de son stage de master 2 à l'Université de Montpellier et modifié par Ambre Petit dans le cadre de son stage de master 1 à l'Université de Montpellier #########################
###########################################################################################################################################
#######          Les versions des outils ici utilisés sont : -python v3.6 et ultérieur                                      ###################
#######                                                      -busco v5.2.2                                              ###################
#######                                                      -augustus v3.5                                             ###################
#######                                                      -bedtools v2.30.0                                          ###################
#######                                                      -Biopython                                                 ###################
#######							     -HMMER v3.3.2                                              ###################
#######							     -ORA v2.0                                                  ###################
#######							     -GRA                                                       ###################
###########################################################################################################################################
###########################################################################################################################################

import yaml
import Bio.SeqIO as SeqIO     #Si l'import ne fonctionne pas, essayer "from Bio import SeqIO"
import os
import psycopg2

# Lecture du fichier de configuration YAML
with open("config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

#Configuration des variables d'environnement 
os.environ["AUGUSTUS_CONFIG_PATH"] = config["augustus_config_path"]
os.environ["PYTHONPATH"] = config["pythonpath"]

# Récupération des paramètres du fichier de configuration
list_species = config["list_species"]
list_DB = config["list_DB"]
genome_dir = config["genome_dir"]
output_dir = config["output_dir"]
database = config["dbname"]
username = config["user"]
hostname = config["host"]
portnumber = config["port"]
extractionline_dir = config["extractionline_dir"]
multifastaToBd_dir = config ["multifastaToBd_dir"]


## Cette règle contient en entrée le dernier fichier généré par le pipeline. Elle permet d'automatiser le lancement de ce dernier.
rule all:
        input:
            expand(output_dir+"/{espece}/complete_{espece}_{DB}.txt" , espece=list_species, DB=list_DB)


## Cette règle permet de changer l’extension du fichier .fna en .fasta pour les gènes provenant du NCBI. 
rule extension:
        input :
                genome_dir+"/{espece}_{DB}.fna"
        output:
                genome_dir+"/{espece}_{DB}.fasta"
        shell :
                "mv {input} {output}"
                

## Cette règle permet d'utiliser l'outil augustus. C'est un outil ab initio qui génère un fichier au format gff contenant les coordonées génomiques des gènes potentiels contenu dans un génome au format fasta.
rule augustus :
        input : 
                genome_dir+"/{espece}_{DB}.fasta"
        output:
                output_dir+"/{espece}/augustus/{espece}_{DB}.gff"
        params: 
                species = config["augustus_species"],
                protein = config["augustus_protein"],
                introns = config["augustus_introns"],
                start = config["augustus_start"],
                stop = config["augustus_stop"],
                cds = config["augustus_cds"],
                exonnames = config["augustus_exonnames"]
        shell : 
                "augustus --species={params.species} --protein={params.protein} --introns={params.introns} --start={params.start} --stop={params.stop} --cds={params.cds} --exonnames={params.exonnames} {input} > {output}"


## Cette règle récupère le fichier généré par la règle précedente et applique un filtre permettant de récupérer seulement les lignes où il est écrit "gene". Cela permet d'éviter les doublons du fichier sur lequel il y a écrit le gène et le transcrit. 
rule augustus_only_gene :  
        input :
                output_dir+"/{espece}/augustus/{espece}_{DB}.gff"
        output:
                output_dir+"/{espece}/augustus/{espece}_{DB}_genes.gff"
        shell : 
                f"python3 {extractionline_dir}/extraction_line.py {{input}} {{output}}"


## La fonction getfasta de bedtools permet de croiser des fichiers au format fasta et gff afin de pouvoir extraire les portions génomiques correspondantes au format fasta. Cette étape est nécessaire en vue de l'utilisation du prochain outil. Elle prend en entrée le fichier fasta du début ainsi que le dernier fichier gff obtenue avec augustus. 
rule bedtools : 
        input: 
                fasta=genome_dir+"/{espece}_{DB}.fasta",
                gff=output_dir+"/{espece}/augustus/{espece}_{DB}_genes.gff"
        output:
                output_dir+"/{espece}/bedtools/{espece}_{DB}_OR.fasta"
        shell: 
                "bedtools getfasta -fi {input.fasta} -bed {input.gff} -fo {output} " 


## La règle split permet de séparer les gènes de taille inférieure et supérieure à une valeur seuil en écrivant dans le fichier lower_length.fasta si inférieur et superior_length.fasta si supérieur.
rule split: 
        input:
                output_dir+"/{espece}/bedtools/{espece}_{DB}_OR.fasta"
        output: 
                output_dir+"/{espece}/bedtools/{espece}_{DB}_OR_lower_length.fasta",
                output_dir+"/{espece}/bedtools/{espece}_{DB}_OR_superior_length.fasta"
        params: 
                max_length = config["split_max_length"]
        run: 
                with open(input[0],"r") as fasta_file :
                        with open (output[0],"w") as sortie:
                                with open (output[1],"w") as output_alternative:
                                        for record in SeqIO.parse(fasta_file,"fasta"): 
                                                if len(record.seq)<= int(params.max_length):
                                                        sortie.write('>'+ str(record.id)+'\n'+ str(record.seq)+ '\n')
                                                else:
                                                        output_alternative.write('>'+ str(record.id)+'\n'+ str(record.seq)+ '\n') 


## ORA est un outil permettant d'identifier les gènes étant des gènes olfactifs. Il utilise en entrée un fichier fasta. 
rule ORA:
        input:
                output_dir+"/{espece}/bedtools/{espece}_{DB}_OR_lower_length.fasta" 
        output:
                output_dir+"/{espece}/ORA/{espece}_{DB}_OR_list.fasta"
        shell : 
                "or.pl --sequence={input} > {output}"


## Cette règle compte le nombre de gène dans le fichier fasta avant analyse, dans le fichier en sortie de ORA ainsi que le nombre de gène qui n'ont pas pu etre analysé par ORA car il possédait une taille supérieure à la limite d'ORA et l'inscrit dans un fichier texte situé à l'emplacement du dossier de l'espèce analysé. 
rule number_of_Gene : 
        input : 
               before=output_dir+"/{espece}/bedtools/{espece}_{DB}_OR.fasta",
               after=output_dir+"/{espece}/ORA/{espece}_{DB}_OR_list.fasta",
               lower=output_dir+"/{espece}/bedtools/{espece}_{DB}_OR_lower_length.fasta",
               upper=output_dir+"/{espece}/bedtools/{espece}_{DB}_OR_superior_length.fasta"
        output:
               output_dir+"/{espece}/{espece}_{DB}_NumberOfGene.txt"
        params: 
                max_length = config["split_max_length"]
        shell :
               """
               total_genes=$(grep -c ">" {input.before} || echo 0)
               or_genes=$(grep -c ">" {input.after} || echo 0)
               superior_genes=$(grep -c ">" {input.upper} || echo 0)
               percentage=$(echo "100 * ($superior_genes / ($total_genes + ($total_genes == 0)))" | bc)

               echo 'Number of all genes before analysis:' > {output} && echo $total_genes >> {output}
               echo 'Number of OR genes after ORA analysis:' >> {output} && echo $or_genes >> {output}
               echo 'Number of genes not analyzed by ORA because of their length greater than {params.max_length} bp:' >> {output} && echo $superior_genes >> {output}
               echo 'Percentage of genes with a length greater than {params.max_length} bp:' >> {output} && echo $percentage >> {output}
               """


## Ajout de l'expérience à la base de données 

connection = psycopg2.connect(dbname=database, user=username, host=hostname, port=portnumber) #connexion à la base de données
new_experience="INSERT INTO experience (\"Pipeline\") VALUES ('Pipeline EXTASOR')"  #ajout de l'experience dans la table experience et génération d'un identifant d'experience
cur=connection.cursor()
cur.execute(new_experience)
connection.commit()

# Récupération de l'identifiant généré 
recup_experience="SELECT max(\"ID\") from experience"
cur.execute(recup_experience)
ID_experience=cur.fetchone()


## Cette règle permet de lancer le script de remplissage de la base de données. Elle émet un fichier temporaire en sortie afin de s'assurer dans "rule all" que la règle ait bien été exécutée.
rule complete_DB :
        input:  
                output_dir+"/{espece}/ORA/{espece}_{DB}_OR_list.fasta"
        output:
                temp(output_dir+"/{espece}/complete_{espece}_{DB}.txt")
        shell :
                f"""
                python3 {multifastaToBd_dir}/multifasta_to_bd.py {{wildcards.espece}} {{wildcards.DB}} {{input}} {{ID_experience}}
                touch {{output}}
                """


## Cette règle permet de supprimer le fichier temporaire créer précedemment. 
rule cleanup:
        input:
                output_dir+"/{espece}/complete_{espece}_{DB}.txt"
        shell:
                "rm {input}"

