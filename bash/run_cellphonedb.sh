me=$(whoami)
source /applications/cpdb-venv/bin/activate

## Use cellphonedb
cellphonedb method statistical_analysis --iterations=1000 --threads=28 \
    --counts-data hgnc_symbol \
    --project-name aa_healthy \
    aa_healthy_meta.txt aa_healthy_counts.txt
