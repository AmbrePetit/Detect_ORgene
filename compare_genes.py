### Comparaison des gènes dans deux fichiers multifasta : celui en sortie du lancement manuel séparé des différents logiciels t celui en sortie du snakemake. 
# file1 : fichier provenant du lancement manuel
# file2 : fichier provenant du snakemake 

# Chemins des fichiers 
file_pairs = [ 
('/home/apetit/Bureau/TestSepare/Pontoporia_blainvillei_NCBI_OR_list.fasta', '/home/apetit/Bureau/TestSepare/Pontoporia_blainvillei_NCBI_NumberOfGene.txt', '/home/apetit/Bureau/Donnees/Pontoporia_blainvillei/ORA/Pontoporia_blainvillei_NCBI_OR_list.fasta', '/home/apetit/Bureau/Donnees/Pontoporia_blainvillei/Pontoporia_blainvillei_NCBI_NumberOfGene.txt','/home/apetit/Bureau/Comparaison/Pontoporia_blainvillei.txt')]

#Boucle sur les paires de fichiers
for file_pair in file_pairs:
	file1_path, file1_numberOfGene, file2_path, file2_numberOfGene, output_file = file_pair


def extract_headers(file_path):
    headers = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                gene_name = line[1:].strip()
                headers.append(gene_name)
    return headers

# Extraire les en-tete des gènes des deux fichiers
headers_file1 = extract_headers(file1_path)
headers_file2 = extract_headers(file2_path)

# Trouver les gènes communs aux deux fichiers
common_genes=[gene for gene in headers_file1 if gene in headers_file2]

# Récupération des lignes dans le fichier qui compte le nombre de gène à la suite du pipeline
with open(file1_numberOfGene, 'r') as file1:
	lines = file1.readlines()
line2_numberOfGene = lines[1].strip()

with open(file2_numberOfGene, 'r') as file2:
	lines = file2.readlines()
line4_numberOfGene = lines[3].strip()

#Calcul du pourcentage d'identité entre les deux fichiers
percentage_common_genes = (int(len(common_genes)) / int(max(line4_numberOfGene, line2_numberOfGene)))*100

# Afficher le nombre de gènes communs
with open(output_file, "w") as output:
        output.write("Pourcentage de genes communs entre les deux fichiers : "+ str(percentage_common_genes)+ "%" +"\nNombre de genes communs : "+str(len(common_genes))+"\n\t sachant qu'il y a "+str(line4_numberOfGene)+" genes dans le fichier issus du snakemake et "+str(line2_numberOfGene)+" genes dans celui issus du lancement manuel du pipeline ")
