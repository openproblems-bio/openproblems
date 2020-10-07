snakemake -j $(grep -c processor /proc/cpuinfo) --use-singularity --singularity-args "-B $(pwd):/opt/openproblems" all
