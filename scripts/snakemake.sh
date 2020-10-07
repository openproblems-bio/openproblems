snakemake -j $(grep -c processor /proc/cpuinfo) --use-singularity --singularity-prefix /tmp/.snakemake --singularity-args "-B $(dirname $(pwd)):/opt/openproblems" all
