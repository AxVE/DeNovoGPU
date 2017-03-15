# DeNovoGPU
L'objectif de ce projet est de tester si les GPU récents peuvent supporter la charge mémoire (avec un algorithme adapté)
d'un assemblage DeNovo.
Le but actuel n'est **pas** de créer un nouvel outil.

## Principe simplifié
Le principe de base est simpliste :
1. On récupère les reads, chaque reads est considéré comme un contig.
2. On calcule un score de fusion entre chaque contig dans une matrice (facilement parallélisable).
3. On fusionne les contigs avec un score supèrieur à une valeur minimale en faisant attention au "cross-fusion"
(2 paires fusionnables ayant un même contig).
4. Si au moins 2 contigs ont pu être fusionnés, on retourne à l'étape 2.
