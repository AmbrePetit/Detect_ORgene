import os
import subprocess
import tempfile
import shutil
import pytest

def test_ORA():
    # Chemins d'entrée et de sortie
    input_file = "/Users/ambre/Desktop/Lipotes_vexillifer/bedtools/Lipotes_vexillifer_NCBI_OR_lower_length.fasta"
    output_file = "/Users/ambre/Desktop/Lipotes_vexillifer/ORA/Lipotes_vexillifer_NCBI_OR_list.fasta"

    # Créer un répertoire temporaire pour exécuter le test
    with tempfile.TemporaryDirectory() as tempdir:
        # Copier le fichier d'entrée dans le répertoire temporaire
        temp_input_file = os.path.join(tempdir, os.path.basename(input_file))
        shutil.copy(input_file, temp_input_file)

        # Exécuter la commande Snakemake dans le répertoire temporaire
        try:
            subprocess.check_output(["snakemake", output_file, "--cores 1"], cwd=tempdir, shell=True)
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Erreur lors de l'exécution de Snakemake : {e.output}")

        # Vérifier si le fichier de sortie existe
        assert os.path.isfile(output_file), "Le fichier de sortie n'a pas été généré"

        # Supprimer le fichier de sortie après le test
        os.remove(output_file)

