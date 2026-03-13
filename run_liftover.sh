export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

for f in /data/cteu/2_users/zyu_lab/mkaminen/overlap_snps/*_hg38.bed; do
    out=${f/_hg38.bed/_hg19.bed}
    echo "Lifting $f..."
    /PHShome/mi415/.local/bin/CrossMap bed hg38ToHg19.over.chain.gz "$f" "$out"

done

