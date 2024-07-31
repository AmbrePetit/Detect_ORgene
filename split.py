import sys

# Vérification du nombre d'arguments
if len(sys.argv) != 3:
    print("Usage: python script.py fichier_entree.fasta fichier_sortie.fasta")
    sys.exit(1)

fichier_entree = sys.argv[1]
fichier_sortie = sys.argv[2]

# Ouverture du fichier fasta en lecture
with open(fichier_entree, "r") as fasta_file:

    # Initialisation des variables pour stocker les informations sur la séquence en cours
    sequence = ""
    header = ""
    length = 0

    # Ouverture du fichier de sortie en écriture
    with open(fichier_sortie, "w") as output_file:

        # Boucle de lecture des lignes du fichier fasta
        for line in fasta_file:

            # Si la ligne est l'en-tête d'une nouvelle séquence
            if line.startswith(">"):

                # Si la séquence précédente est un gène court, on l'écrit sur le fichier de sortie
                if length < 2500 and sequence != "":
                    output_file.write(header + sequence + "\n")

                # On initialise les variables pour la nouvelle séquence
                header = line
                sequence = ""
                length = 0

            # Sinon, on ajoute la ligne à la séquence en cours et on met à jour la longueur de la séquence
            else:
                sequence += line.strip()
                length = len(sequence)

        # À la fin de la boucle, on écrit la dernière séquence si elle est un gène court
        if length < 2500 and sequence != "":
            output_file.write(header + sequence + "\n")
