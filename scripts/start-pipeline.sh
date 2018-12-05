REFERENCES=/home/oskar/01-workspace/01-data/refdata/hg38
docker run --rm -u $UID:1000 -ti -v $REFERENCES:/references -v `pwd`:/data -w /data oskarv/snakemake-germline-tools snakemake -j -T
