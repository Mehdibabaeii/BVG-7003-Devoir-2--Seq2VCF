#!/bin/bash

printf "\nCréation de la structure des répertoires : data, barcodes, refgenome, results, logs\n"

# Création des dossiers principaux
for i in data barcodes refgenome results logs
do
    mkdir -p $i
done

# Création des sous-dossiers dans 'results'
printf "\nAjout des sous-dossiers dans 'results' : alignment, demultiplexed, qc_reports, trimmed, variants\n"
for subdir in alignment demultiplexed qc_reports trimmed variants
do
    mkdir -p results/$subdir
done
