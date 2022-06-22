# STAGE-2022

## NOTE : TOUTES LES OPÉRATIONS SONT EFFECTUÉES SUR LES CLUSTERS CALCUL QUÉBEC

Pour préparer les fichier .bg :

1) Copier les fichiers **smart2_cov.py** et **matrix.sh** dans le dossier avec les fichiers .cov.gz bien renommés.
2) Lancer matrix.sh : le script va créer un dossier ./matrix avec tous les fichiers .bg dedans. Effectuer la commande, montrée en exemple ci-dessous, avec les fichiers résultants.

```bash
module load bedtools/2.30.0
bedtools unionbedg -i G1-1.bg G1-2.bg G2-1.bg G2-2.bg -header -names G1_1 G1_2 G2_1 G2_2 -filler - > MethylMatrix.txt &
```

<details>
  <summary>Exemple complet</summary>
  
### Pipeline du labo + Genpipes
  
  ```bash
module load bedtools/2.30.0
bedtools unionbedg -i CT-M16-1.bg CT-M16-2.bg CT-M16-3.bg CT-M16-4.bg HFD-M16-1.bg HFD-M16-2.bg HFD-M16-3.bg HFD-M16-4.bg CT-F16-1.bg CT-F16-2.bg CT-F16-3.bg HFD-F16-1.bg HFD-F16-2.bg HFD-F16-3.bg CT-M18-1.bg CT-M18-2.bg CT-M18-3.bg CT-M18-4.bg -header -names CT-M16_1 CT-M16_2 CT-M16_3 CT-M16_4 HFD-M16_1 HFD-M16_2 HFD-M16_3 HFD-M16_4 CT-F16_1 CT-F16_2 CT-F16_3 HFD-F16_1 HFD-F16_2 HFD-F16_3 CT-M18_1 CT-M18_2 CT-M18_3 CT-M18-4 -filler - > MethylMatrix.txt &
```

</details>

Lancer SMART2, en supposant que l'environnement virtuel est bien configuré :

```bash
module load StdEnv/2018.3 python/2.7 scipy-stack
source ~/ENV/bin/activate
SMART (Matrice) -t DeNovoDMR -o ./(dossier output) -AD 0.10 -PC 1
```

Dans ce cas-ci, on suppose que l'environnement virtuel se trouve à la racine pour venir l'activer. L'option -AD est la différence de méthylation moyenne absolue entre le groupe de cas et le groupe de contrôle, ajustée à 10%. L'option -PC est la **p-value** auquels les DMR des case-controls sont identifiés.
