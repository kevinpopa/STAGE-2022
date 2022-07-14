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
  
Pipeline du labo + Genpipes
  
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

## Analyse des segments - Taux de méthylation
À partir du fichier 3_MergedSegment_Methylation fourni par SMART2, on veut déterminer la méthylation moyenne de chaque réplicat.

```R
segment_file <- read.delim("C:/Users/Admin/Desktop/3_MergedSegment_Methylation.txt", comment.char="#")
# On veut seulement que la dernière colonne
sample_methylation <- segment_file["Sample_Methyl.CT.M16_1.CT.M16_2.CT.M16_3.CT.M16_4.HFD.M16_1.HFD.M16_2.HFD.M16_3.HFD.M16_4.CT.F16_1.CT.F16_2.CT.F16_3.HFD.F16_1.HFD.F16_2.HFD.F16_3.CT.M18_1.CT.M18_2.CT.M18_3.CT.M18_4"]

# Ensuite, on extrait cette dernière colonne vers un fichier texte, on re-importe les données sous forme de fichier csv 
write.table(sample_methylation, file = "sample_methylation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
sample_meth <- read.csv("~/Segment Methylation/sample_methylation.txt")

# On s'intéresse à la moyennne de chaque replicat (colonne de simple_meth)
means_samples <- colMeans(as.matrix(sapply(sample_meth, as.numeric)), na.rm = TRUE)
means_samples

 CT.M16_1   CT.M16_2   CT.M16_3   CT.M16_4  HFD.M16_1  HFD.M16_2  HFD.M16_3  HFD.M16_4   CT.F16_1   CT.F16_2   CT.F16_3 
0.18669140 0.15081660 0.13311596 0.11500765 0.10781484 0.09872396 0.15331439 0.13950453 0.08410354 0.04107823 0.03127384 
 HFD.F16_1  HFD.F16_2  HFD.F16_3   CT.M18_1   CT.M18_2   CT.M18_3   CT.M18_4 
0.04584711 0.04148280 0.03915652 0.03741029 0.02782513 0.05099672 0.04286993 

# On effectue un t-test pour voir si la différence des moyennes de méthylation des deux groupes est significativement différente. 
# Ici, en exemple, on vérifie si la moyenne de méthylation des témoins mâles (1 à 4) est différente des traitement HFD mâles (5 à 8). 
# Comme vu ci-dessus, une p-value de 0.3235 est non significative. 
t.test(means_samples[1:4], means_samples[5:8], alternative = "two.sided", var.equal = FALSE)

	Welch Two Sample t-test

data:  means_samples[1:4] and means_samples[5:8]
t = 1.0781, df = 5.8352, p-value = 0.3235
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.02772006  0.07085701
sample estimates:
mean of x mean of y 
0.1464079 0.1248394 
```

<details>
  <summary>Les deux autres t-test</summary>
  
Témoins femelles (9 à 11) vs. Traitement HFD femelles (12 à 14) & Témoins mâles (1 à 4) vs. Témoins femelles (9 à 11)
  
  ```R
> t.test(means_samples[9:11], means_samples[12:14], alternative = "two.sided", var.equal = FALSE)

	Welch Two Sample t-test

data:  means_samples[9:11] and means_samples[12:14]
t = 0.61127, df = 2.0584, p-value = 0.6017
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.05844882  0.07842827
sample estimates:
 mean of x  mean of y 
0.05215187 0.04216214 

> t.test(means_samples[1:4], means_samples[9:11], alternative = "two.sided", var.equal = FALSE)

	Welch Two Sample t-test

data:  means_samples[1:4] and means_samples[9:11]
t = 4.2281, df = 4.6726, p-value = 0.009584
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.03572113 0.15279093
sample estimates:
 mean of x  mean of y 
0.14640790 0.05215187
```

</details>

```R
#On peut storer les moyennes des groupes dans des variables 
CTM_mean <- mean(means_samples[1:4])
HFDM_mean <- mean(means_samples[5:8])
CTF_mean <- mean(means_samples[9:11])
HFDF_mean <- mean(means_samples[12:14])

#On produit un barplot. Pour arriver à cela, on crée une table ayant pour colonne chaque échantillon. 
tab_means <- matrix(c(CTM_mean, CTF_mean, HFDM_mean, HFDF_mean), ncol = 4, byrow = TRUE)
rownames(tab_means) <- c("mean")
colnames(tab_means) <- c("Témoin mâle (n=4)", "Témoin femelle (n=3)", "HFD mâle (n=4)", "HFD femelle (n=3)")
tab_means <- as.table(tab_means)
barplot(tab_means, ylim = c(0,1))

```

## Analyse des DMR - Homer
Avec les fichiers pour chaque cas d'hypométhylation et hyperméthylation, on va vouloir annoter les différents DMR sur le génome du rat. 

À partir du fichier script en R (script_prepareHomer.R : code ci-dessous), préparer les fichiers pour l'annotation Homer. 

```bash
module load r/4.1.2
Rscript script_prepareHomer.R (fichier input) (fichier output : optionnel)
```

```R
args = commandArgs(trailingOnly = TRUE)

# Check that at least the input file is specified.
if (length(args)==0) {
  stop("The script input should be Rscript script_prepareHomer.R (input file : required) (output file name : optional)", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(gsub('.{4}$', '', args[1]), "HOMER.txt", sep = ".")
}
df <- read.delim(args[1])

df_1 <- df[, c("RegionChrome", "RegionStart", "RegionEnd")]
df_1$peak_name <- paste(df_1$RegionChrome, df_1$RegionStart, df_1$RegionEnd, sep = "_")
df_1$NV <- "NV"
df_1$Strand <- "*"

write.table(df_1, args[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
