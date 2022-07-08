# STAGE-2022

## NOTE : TOUTES LES OPÉRATIONS SONT EFFECTUÉES SUR LES CLUSTERS CALCUL QUÉBEC

Pour préparer les fichier .bg :

1) Copier les fichiers **smart2_cov.py** et **matrix.sh** dans le dossier avec les fichiers .cov.gz bien renommés.
2) Lancer matrix.sh : le script va créer un dossier ./matrix avec tous les fichiers .bg dedans. Effectuer la commande, montrée en exemple ci-dessous, avec les fichiers résultants.

_Note : Bedtools reconnait les réplicats de groupes que s'ils sont nommés avec des barres en bas._ 

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

Ensuite, à partir du fichier "7_DMR_Case_control.txt", qui est un tableau ayant tous les DMR identifiés avec les paramètres voulus, séparer les différents DMR pour chaque duo de cas-contrôle. Avant de créer tous les fichiers voulus, on veut effectuer une correction pour tests multiples "FDR".

Ce code en R permet d'importer les données et ajouter une colonne FDR au dataframe (df). 

```R
df <- read.delim("C:/Users/Admin/Desktop/7_DMR_Case_control.txt", comment.char="#")
p <- df$TwoSampleTtest_pvalue
qvalue <- p.adjust(p, method = "fdr", n = length(p))
df['FDR'] <- qvalue
```

Ensuite, on veut que seulement les DMR ayant une q-value inférieure ou égale à 0.05. 

```R
dfsig <- df[df$FDR <= 0.05, ]
```

Finalement, on séparera ce dataframe en plusieurs sous-dataframes pour les différentes conditions de cas-contrôle et l'état de la méthylation. Ces étapes pourraient être simplifiées et combinées pour moins de ligne de code, mais voici la version longue. 

```R
dfCTM_HFDM <- dfsig[dfsig$Case_Group == "CT-M16" & dfsig$Control_Group == "HFD-M16", ]
dfCTM_HFDM <- dfsig[dfsig$Case_Group == "HFD-M16" & dfsig$Control_Group == "CT-M16", ]
dfCTF_HFDF <- dfsig[dfsig$Case_Group == "HFD-F16" & dfsig$Control_Group == "CT-F16", ]
dfCTM_CTF <- dfsig[dfsig$Case_Group == "CT-M16" & dfsig$Control_Group == "CT-F16", ]
dfCT18_CTM <- dfsig[dfsig$Case_Group == "CT-M18" & dfsig$Control_Group == "CT-M16", ]
dfCT18_CTF <- dfsig[dfsig$Case_Group == "CT-M18" & dfsig$Control_Group == "CT-F16", ]

CTM_HFDM_Hyper <- dfCTM_HFDM[dfCTM_HFDM$Case_Status == "Hyper", ]
CTM_HFDM_Hypo <- dfCTM_HFDM[dfCTM_HFDM$Case_Status == "Hypo", ]
CTF_HFDF_Hyper <- dfCTF_HFDF[dfCTF_HFDF$Case_Status == "Hyper", ]
CTF_HFDF_Hypo <- dfCTF_HFDF[dfCTF_HFDF$Case_Status == "Hypo", ]
CTM_CTF_Hyper <- dfCTM_CTF[dfCTM_CTF$Case_Status == "Hyper", ]
CTM_CTF_Hypo <- dfCTM_CTF[dfCTM_CTF$Case_Status == "Hypo", ]
CT18_CTM_Hyper <- dfCT18_CTM[dfCT18_CTM$Case_Status == "Hyper", ]
CT18_CTM_Hypo <- dfCT18_CTM[dfCT18_CTM$Case_Status == "Hypo", ]
CT18_CTF_Hyper <- dfCT18_CTF[dfCT18_CTF$Case_Status == "Hyper", ]
CT18_CTF_Hypo <- dfCT18_CTF[dfCT18_CTF$Case_Status == "Hypo", ]
```

Utiliser pour chacun des dataframes résultants le code suivant pour l'exporter en fichier texte. Ceci produit un fichier dans le répertoire de travail.

```R
write.table(CTM_HFDM_Hyper,"CTM_HFDM_Hyper.txt", sep= "\t", quote = FALSE)
```
