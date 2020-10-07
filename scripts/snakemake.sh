snakemake -j 1 --use-singularity --singularity-prefix /tmp/.snakemake --singularity-args "-B $(dirname $(pwd)):/opt/openproblems" all
