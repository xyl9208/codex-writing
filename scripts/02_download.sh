#!/bin/bash
# Download processed expression data + platform annotations from GEO FTP
set -u
cd "$(dirname "$0")/../data/geo"
dl() { url=$1; out=$(basename "$url"); if [ -s "$out" ]; then echo "exists $out"; else curl -sS -L --retry 4 --retry-delay 3 -o "$out" "$url" && echo "ok $out $(du -h "$out" | cut -f1)" || echo "FAIL $url"; fi; }
G=https://ftp.ncbi.nlm.nih.gov/geo/series
dl $G/GSE26nnn/GSE26155/matrix/GSE26155_series_matrix.txt.gz
dl $G/GSE140nnn/GSE140947/suppl/GSE140947_bulkrna_24_samp_genecounts.csv.gz
dl $G/GSE235nnn/GSE235161/suppl/GSE235161_expr_fpkm_42sample.csv.gz
dl $G/GSE202nnn/GSE202267/suppl/GSE202267_GEO_TAA_AAA_readcounts.csv.gz
dl $G/GSE52nnn/GSE52093/matrix/GSE52093_series_matrix.txt.gz
dl $G/GSE52nnn/GSE52093/suppl/GSE52093_non-normalized.txt.gz
dl $G/GSE153nnn/GSE153434/suppl/GSE153434_all.counts.txt.gz
dl $G/GSE98nnn/GSE98770/matrix/GSE98770-GPL14550_series_matrix.txt.gz
dl $G/GSE190nnn/GSE190635/matrix/GSE190635_series_matrix.txt.gz
dl $G/GSE219nnn/GSE219204/matrix/GSE219204_series_matrix.txt.gz
dl $G/GSE147nnn/GSE147026/suppl/GSE147026_mRNA-ADVSCK.All.anno.txt.gz
dl $G/GSE294nnn/GSE294606/suppl/GSE294606_proceeded_data.txt.gz
dl $G/GSE267nnn/GSE267434/suppl/GSE267434_transcript_sample_FPKM.txt.gz
dl $G/GSE318nnn/GSE318877/suppl/GSE318877_RAW.tar
dl $G/GSE314nnn/GSE314998/suppl/GSE314998_Group_one.txt.gz
dl $G/GSE314nnn/GSE314998/suppl/GSE314998_Group_two.txt.gz
for s in GSE318877 GSE155468 GSE213740 GSE189795 GSE219204; do stub=${s:0:$((${#s}-3))}nnn; curl -sS -L -o ${s}_filelist.txt $G/$stub/$s/suppl/filelist.txt && echo "ok ${s}_filelist.txt"; done
P=https://ftp.ncbi.nlm.nih.gov/geo/platforms
dl $P/GPL5nnn/GPL5175/annot/GPL5175.annot.gz
dl $P/GPL10nnn/GPL10558/annot/GPL10558.annot.gz
dl $P/GPL14nnn/GPL14550/annot/GPL14550.annot.gz
dl $P/GPL570/GPL570/annot/GPL570.annot.gz
dl $P/GPL24nnn/GPL24539/annot/GPL24539.annot.gz
ls -la
