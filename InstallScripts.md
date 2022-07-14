Pour avoir accès aux fichiers de scripts partout dans Linux, on va les déplacer dans un endroit dans lequel Linux reconnait automatiquement les chemin vers ces fichiers. 
Si non créé, créer un dossier bin à la racine.

```bash
mkdir ~/bin
```

Rajouter les fichiers matrix.sh et smart2_cov.py dans ~/bin.

Donnez les permissions de lecture et exécution au fichier pour tout le monde (si le fichier n'est pas final, rajouter l'écriture).

```bash
chmod ugo+rx matrix.sh
```

Maintenant, il est possible de faire matrix.sh n'importe où sans spécifier son chemin complet à chaque fois. Faire attention de seulement lancer le script dans le bon dossier avec les cov.gz bien renommés. Il est aussi conseillé de lancer le script en arrière-plan comme suit :

```bash
matrix.sh &
```
