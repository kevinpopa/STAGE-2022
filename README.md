# STAGE-2022

Pour préparer les fichier .bg :

1) Copier les fichiers **smart2_cov.py** et **matrix.sh** dans le dossier avec les fichiers .cov.gz bien renommés.
2) Lancer matrix.sh : le script va créer un dossier ./matrix avec tous les fichiers .bg dedans. Effectuer la commande ci-dessous avec les fichiers résultants.

Voici un exemple :

```bash
module load bedtools/2.30.0
bedtools unionbedg -i G1-1.bg G1-2.bg G2-1.bg G2-2.bg -header -names G1-1 G1-2 G2-1 G2-2 -filler - > MethylMatrix.txt &
```

Lancer SMART2, en supposant que l'environnement virtuel est bien configuré :

```bash
module load StdEnv/2018.3 python/2.7 scipy-stack
source ~/ENV/bin/activate
SMART (Matrice) -t DeNovoDMR -o ./(dossier output) -AD 0.10
```
