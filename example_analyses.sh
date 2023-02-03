!/bin/bash

genomes="\
Cinnamomum_micranthum_f_kanehirae  Laurales     1  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/546/025/GCA_003546025.1_ASBRC_Ckan_1.0/GCA_003546025.1_ASBRC_Ckan_1.0_genomic.fna.gz
Pyrus_betulifolia                  Rosales      0  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/844/245/GCA_007844245.1_ASM784424v1/GCA_007844245.1_ASM784424v1_genomic.fna.gz
Setaria_italica                    Poales       1  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/263/155/GCA_000263155.2_Setaria_italica_v2.0/GCA_000263155.2_Setaria_italica_v2.0_genomic.fna.gz
Triticum_aestivum                  Poales       0  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/294/505/GCA_018294505.1_IWGSC_CS_RefSeq_v2.1/GCA_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz
Aegilops_tauschii_ssp_strangulata  Poales       0  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/575/655/GCA_002575655.2_Aet_v5.0/GCA_002575655.2_Aet_v5.0_genomic.fna.gz
Coffea_arabica                     Gentianales  0  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/713/225/GCA_003713225.1_Cara_1.0/GCA_003713225.1_Cara_1.0_genomic.fna.gz
Coffea_canephora                   Gentianales  1  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/059/795/GCA_900059795.1_AUK_PRJEB4211_v1/GCA_900059795.1_AUK_PRJEB4211_v1_genomic.fna.gz
Coffea_eugenioides                 Gentianales  0  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/713/205/GCA_003713205.1_Ceug_1.0/GCA_003713205.1_Ceug_1.0_genomic.fna.gz
Artemisia_annua                    Asterales    1  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/112/345/GCA_003112345.1_ASM311234v1/GCA_003112345.1_ASM311234v1_genomic.fna.gz
Cynara_cardunculus_var_scolymus    Asterales    1  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/531/365/GCA_001531365.2_CcrdV1.1/GCA_001531365.2_CcrdV1.1_genomic.fna.gz"

# download and prepare genomic data
while IFS=$'\n' read genome; do
    arr=($genome)
    taxon=${arr[0]}
    annotated=${arr[2]}
    url=${arr[3]}
    
    # download genome
    wget --read-timeout 30 "$url"
    
    # download annotation
    if [ "$annotated" -eq 1 ]; then
        url_ann="$(sed 's/fna\.gz$/gff.gz/' <<< $url)"
        wget --read-timeout 30 "$url_ann"
    fi
    
    # verify checksums
    wget "$(dirname $url)/md5checksums.txt"
    echo -e "\ncheck files:" 
    md5sum --ignore-missing -c md5checksums.txt
    if [ $? -ne 0 ]; then
        echo "Error: corrupted download!"
        exit
    fi
    echo ""
    rm -f md5checksums.txt
    
    # extract and rename files
    gunzip  "$(basename $url)"
    mv "$(basename $url | sed 's/\.gz$//')" "${taxon}.fna"
    if [ $annotated -eq 1 ]; then
        url_ann="$(sed 's/fna\.gz$/gff.gz/' <<< $url)"
        gunzip  "$(basename $url_ann)"
        mv "$(basename $url_ann | sed 's/\.gz$//')" "${taxon}.gff"
    fi
done <<< "$genomes"


while IFS=$'\n' read genome; do
    arr=($genome)
    echo ${arr[0]}, ${arr[1]}, ${arr[2]}, ${arr[3]}
done <<< "$genomes"

# download and prepare target sequences
wget https://github.com/mossmatters/Angiosperms353/raw/master/Angiosperms353_targetSequences.fasta
wget https://raw.githubusercontent.com/Smithsonian/Compositae-COS-workflow/master/COS_sunf_lett_saff_all.fasta

# problem: lines of different lengths -> normalize
./docker.sh fasta_formatter -i Angiosperms353_targetSequences.fasta -o tmp.fasta -w 100
mv tmp.fasta Angiosperms353_targetSequences.fasta

# download information about sequence IDs
wget https://github.com/mossmatters/Angiosperms353/raw/master/support_files/fine_annotations_MJ.txt
cp fine_annotations_MJ.txt abbreviations.txt

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

# Analyses based on Angiosperms353
##################################

# 1) using whole genome assemblies:

while IFS=$'\n' read genome; do
    arr=($genome)
    taxon=${arr[0]}
    order=${arr[1]}
    annotated=${arr[2]}
    
    rm "${order}.fasta" "${order}.fasta.fai" 2> /dev/null

    gawk -v "o=$order" '
    ARGIND==1 && NR>1 {
        order[$1] = $2
    }
    ARGIND==2 {
        if ($0 ~ /^>/) {
            ID = $0
            gsub(/^>/, "", ID)
            gsub(/-.*/, "", ID)
            if (order[ID] == o) {
                p = 1
            }
            else {
                p = 0
            }
        }
        if (p == 1) {
            print
        }
    }
    ' abbreviations.txt Angiosperms353_targetSequences.fasta > "${order}.fasta"
    
    if [ $annotated -eq 1 ]; then
        time ./docker.sh ./run.py "${order}.fasta" "${taxon}.fna" -a "${taxon}.gff"
    else
        time ./docker.sh ./run.py "${order}.fasta" "${taxon}.fna"
    fi
done <<< "$genomes"


# 2) using single subgenomes from polyploids:

# extract subgnomes (incl. unassigned scaffolds) from allo-tetraploid coffee
fasta_formatter -i Coffea_arabica.fna | awk '/^>/ && !/chromosome [0-9]+e/ {p=1; print} !/^>/ {if (p==1) {print} p=0}' | fasta_formatter -w 100 > Coffea_arabica_c.fna
fasta_formatter -i Coffea_arabica.fna | awk '/^>/ && !/chromosome [0-9]+c/ {p=1; print} !/^>/ {if (p==1) {print} p=0}' | fasta_formatter -w 100 > Coffea_arabica_e.fna

# extract subgnomes (incl. unassigned scaffolds) from allo-hexaploid wheat
fasta_formatter -i Triticum_aestivum.fna | awk '/^>/ && !/chromosome [1-7][BD]/ && !/scaffold[1-7][BD]/ {p=1; print} !/^>/ {if (p==1) {print} p=0}' | fasta_formatter -w 100 > Triticum_aestivum_A.fna
fasta_formatter -i Triticum_aestivum.fna | awk '/^>/ && !/chromosome [1-7][AD]/ && !/scaffold[1-7][AD]/ {p=1; print} !/^>/ {if (p==1) {print} p=0}' | fasta_formatter -w 100 > Triticum_aestivum_B.fna
fasta_formatter -i Triticum_aestivum.fna | awk '/^>/ && !/chromosome [1-7][AB]/ && !/scaffold[1-7][AB]/ {p=1; print} !/^>/ {if (p==1) {print} p=0}' | fasta_formatter -w 100 > Triticum_aestivum_D.fna

for subgenome in c e; do
    time ./docker.sh ./run.py Gentianales.fasta "Coffea_arabica_${subgenome}.fna"
done

for subgenome in A B D; do
    time ./docker.sh ./run.py Poales.fasta "Triticum_aestivum_${subgenome}.fna"
done



# Analyses based on CompCOS
###########################

grep -E ">.{9}sunf" COS_sunf_lett_saff_all.fasta -A1 | grep -v "^--$" > compcos_sunf.fasta
grep -E ">.{9}saff" COS_sunf_lett_saff_all.fasta -A1 | grep -v "^--$" > compcos_saff.fasta

for e in 1e-10 1e-5 1 10; do
    for m in 10 15 25; do
        time ./docker.sh ./run.py -e $e -m $m compcos_sunf.fasta Artemisia_annua.fna -a Artemisia_annua.gff
        time ./docker.sh ./run.py -e $e -m $m compcos_saff.fasta Cynara_cardunculus_var_scolymus.fna -a Cynara_cardunculus_var_scolymus.gff
    done
done

# summarize results
###################

find . -name summary.txt | awk -F"/" '
BEGIN {
    OFS = "\t"
}
{
    summary = $0
    targets = $2
    genome = $3
    e = $4
    m = $5

    sub(/^e_thresh_/, "", e)
    sub(/^mc_thresh_/, "", m)
    
    "grep -cv \"^$\" "summary | getline n_loci
    
    if (genome_size[genome] == 0) {
        "grep -v \"^>\" "genome".fna | tr -d [:space:] | wc -c" | getline genome_size[genome]
    }
    
    "awk '"'"'{i+=$2} END {print i}'"'"' "summary | getline total_length
    
    summary_groups = summary
    sub(/summary\.txt$/, "summary_groupwise.txt", summary_groups)
    grouped = getline < summary_groups < 0 ? 0 : 1
    
    if (grouped == 1) {
        "grep -cv \"^$\" "summary_groups | getline n_groups
        "awk '"'"'{i+=$2} END {print i}'"'"' "summary_groups | getline total_group_length
    }
    else {
        n_groups = 0
        total_group_length = 0
    }
    
    avg_locus_length = n_loci > 0 ? total_length/n_loci : "-"
    avg_group_length = n_groups > 0 ? total_group_length/n_groups : "-"
    
    print targets, genome, genome_size[genome], e, m, n_loci, avg_locus_length, n_groups, avg_group_length
}
' | sort -k1,1 -k2,2 -k5,5g -k4,4g | column -t > results.txt
