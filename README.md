# DeNovoGPU
L'objectif de ce projet est de tester si les GPU récents peuvent supporter la charge mémoire (avec un algorithme adapté)
d'un assemblage DeNovo.
Le but actuel n'est **pas** de créer un nouvel outil.
Il existe déjà des outils (généralement basés sur des graphes de De-Bruijn) exploitant à certains moments la technologie GPU.
Une liste non-exhaustive :
- Velvet
- SOAP
- Abyss
- Ray

## Pré-requis
- **g++:** g++ (GCC) 6.3.1 20170109 (support de c++14)
- **OpenCL:** 2.1
- **OpenCL C++ wrapper API:** 1.2

## Principe simplifié
Le principe de base est simpliste :
1. On récupère les reads, chaque reads est considéré comme un contig.
2. On calcule un score de fusion entre chaque contig dans une matrice (facilement parallélisable).
3. On fusionne les contigs avec un score supèrieur à une valeur minimale en faisant attention au "cross-fusion"
(2 paires fusionnables ayant un même contig).
4. Si au moins 2 contigs ont pu être fusionnés, on retourne à l'étape 2.

## Mesures
Les mesures rapides de temps et d'empreinte mémoire (CPU et RAM) sont réalisées avec l'outil time (version GNU time 1.7).
Attention, la plupart des configurations pointent la commande _time_ sur un vieille version.
Vous pouvez essayer _/usr/bin/time_.

Un exemple d'utilisation :

`/usr/bin/time -f "\nTime:\t%E\nMem:\t%M KB" ./Bin/denovoGPU -f Data/reads_test.fasta -t 6`

## Liste de tâches
- [x] Test en multithreading
- [ ] Implémentation OpenCL et test
- [ ] Nettoyage du code et commentaires
- [ ] Output des contigs finaux
- [ ] Retravaille des contigs finaux pour avoir plus de précisions (en mappant les reads du contigs sur le contigs puis travail par profondeur).
- [ ] Stocker et manipuler les nucléotides des contigs sur 4 bits (soit 2 par char) plutôt que 1 par char.
Cela permet de réduire la charge mémoire sans perte d'information (IUPAC code) : chaque bit d'un 4 bits représente la possibilité d'avoir tel nucléotide à cet position.
Par exemple :
	- 1000: A
	- 0100: G
	- 0010: C
	- 0001: T
	- 1111: X
	- 1011: H
	- ...
- [ ] Multiple formats fichiers entrants (détection automatique)

