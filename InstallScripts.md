Pour avoir accès aux fichiers de scripts partout dans Linux, on va les déplacer dans un endroit dans lequel Linux reconnait automatiquement les chemin vers ces fichiers. 
Si non créé, créer un dossier bin à la racine.

```bash
mkdir ~/bin
```

Rajouter les fichiers dans le dossier scripts de ce Github (STAGE-2022/scripts/) dans ~/bin. Ainsi, le fichier bin devrait ressembler à cela :

```bash
.
├── bin
│   ├── matrix.sh
│   ├── script_prepareHomer.R
│   └── smart2_cov.py
```

Donnez les permissions de lecture et exécution aux fichier pour tout le monde (si le fichier n'est pas final, rajouter l'écriture).  À faire surtout pour les fichiers script shell .sh, pour les rendre exécutables sans erreurs de permissions. 

```bash
chmod ugo+rx ~/bin/matrix.sh
```

Maintenant, il est possible de faire matrix.sh n'importe où sans spécifier son chemin complet à chaque fois. Faire attention de seulement lancer le script dans le bon dossier avec les cov.gz bien renommés. Il est aussi conseillé de lancer le script en arrière-plan comme suit :

```bash
matrix.sh &
```
