DROP TYPE IF EXISTS nvx CASCADE;
DROP TYPE IF EXISTS cptm CASCADE;
DROP TYPE IF EXISTS etats CASCADE;

DROP TABLE IF EXISTS gene CASCADE;
DROP TABLE IF EXISTS experience CASCADE;
DROP TABLE IF EXISTS assemblie CASCADE;
DROP TABLE IF EXISTS link CASCADE;
DROP TABLE IF EXISTS organism CASCADE;

--Creation des types utilisé dans les tables

CREATE TYPE nvx as enum ('Scaffold', 'complet'); 
CREATE TYPE etats as enum ('pseudogène','fonctionnel');


--Creation des tables et relations

CREATE TABLE organism(
	"Genre" varchar(100),
	"Espèce" varchar(100),
	"Micro ordre" varchar(100),
	"Ordre" varchar(100),
	"Habitat" varchar(100),
	"Alimentation" varchar(100),
	"Comportement social" varchar(100),
	"Région" varchar(150),
	CONSTRAINT 
		pk_organisme PRIMARY KEY ("Espèce","Genre")
);


CREATE TABLE assemblie(
	"ID" varchar(200) PRIMARY KEY,
	"Genre" varchar (100), 
	"Espèce" varchar (100),
	"Base de donnee" varchar(100),
	"Date de publication" date,
	"Niveau d'assemblage" nvx,
	"Score_busco" varchar(120),
	CONSTRAINT
		fk_esp FOREIGN KEY ("Espèce","Genre") REFERENCES organism ("Espèce","Genre")
);


CREATE TABLE experience(
	"ID" SERIAL PRIMARY KEY,
	"Pipeline" varchar(200),
	"Paramétres" text
);

CREATE TABLE gene(
	"ID" int PRIMARY KEY,
	"Nom" text,
	"Superfamille" varchar(100),
	"Famille" varchar(10),
	"Etat" etats,
	"Début" int,
	"Fin" int,
	"Séquence" text,
	"Référence" int
);

CREATE TABLE link(
	"ID_assemblie" varchar(200),
	"ID_experience" int,
	"ID_gene" int,
CONSTRAINT 
	pk_link PRIMARY KEY ( "ID_assemblie", "ID_experience","ID_gene"),
CONSTRAINT
	fk_link_assemblie FOREIGN KEY ("ID_assemblie") REFERENCES assemblie("ID"),
CONSTRAINT 
	fk_link_experience FOREIGN KEY ("ID_experience") REFERENCES experience("ID"),
CONSTRAINT 
	fk_link_gene FOREIGN KEY ("ID_gene") REFERENCES gene("ID")
);


