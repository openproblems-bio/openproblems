snakemake -j $(grep -c processor /proc/cpuinfo) --use-singularity --singularity-args "-B $(dirname $(pwd)):/opt/openproblems" all
