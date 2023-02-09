#!/bin/bash

# download and prepare target sequences
wget https://github.com/mossmatters/Angiosperms353/raw/master/Angiosperms353_targetSequences.fasta

# problem: lines of different lengths -> normalize
fasta_formatter -i Angiosperms353_targetSequences.fasta -o tmp.fasta -w 100
mv tmp.fasta Angiosperms353_targetSequences.fasta

# download information about sequence IDs
wget https://github.com/mossmatters/Angiosperms353/raw/master/support_files/fine_annotations_MJ.txt
tail -n+2 fine_annotations_MJ.txt > abbreviations.txt

# add missing information about some abbreviations
# (species assignment according to https://github.com/GrosseLab/OneKP-gene-family-evo/blob/master/data/phylogeny/species_dictionary.tsv)
echo "
Ambtr  Amborellales  Amborellaceae  Amborella    Amborella_trichopoda
Aquco  Ranunculales  Ranunculaceae  Aquilegia    Aquilegia_coerulea
Arath  Brassicales   Brassicaceae   Arabidopsis  Arabidopsis_thaliana
Carpa  Brassicales   Caricaceae     Carica       Carica_papaya
Eucgr  Myrtales      Myrtaceae      Eucalyptus   Eucalyptus_grandis
Manes  Malpighiales  Euphorbiaceae  Manihot      Manihot_esculenta
Mimgu  Lamiales      Phrymaceae     Mimulus      Mimulus_guttatus
Orysa  Poales        Poaceae        Oryza        Oryza_sativa
Phavu  Fabales       Fabaceae       Phaseolus    Phaseolus_vulgaris
Poptr  Malpighiales  Salicaceae     Populus      Populus_trichocarpa
Prupe  Rosales       Rosaceae       Prunus       Prunus_persica
Solly  Solanales     Solanaceae     Solanum      Solanum_lycopersicum
Sorbi  Poales        Poaceae        Sorghum      Sorghum_bicolor
Theca  Malvales      Malvaceae      Theobroma    Theobroma_cacao
Vitvi  Vitales       Vitaceae       Vitis        Vitis_vinifera
" | sed 's/ \+/\t/g' >> abbreviations.txt

sort -k2,2 -k3,3 -k4,4 -k5,5 abbreviations.txt | grep -v "^$" > abbreviations.tmp
mv abbreviations.tmp abbreviations.txt
