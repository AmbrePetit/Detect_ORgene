import psycopg2

con = psycopg2.connect(    #paramètres de connexion à la base de donnée
    dbname='CeGeC',        #nom de la base de donnée
    user='Ambre',          #nom de l'utilisateur
    host='localhost',      #spécifier l'hôte local ou l'adresse IP
    port='5432')           #numéro du port, par défaut 5432
cur = con.cursor()         #création d'un curseur permettant l'excécution de requêtes SQL
cur.execute("requête SQL") #requête SQL que l'on souhaite effectuer
rows = cur.fetchall()      #récupère tous les résultats de la requête précédente dans une liste

for row in rows:         #affichage des resultats
    print(row)

cur.close()                #fermeture du curseur
con.close()                #fermeture de la connexion à la base de données PostgreSQL



