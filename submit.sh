#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=32G,h_rt=20:00:00,tmem=32G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


FOLDER=submissions/$(date +"%Y%m%d%H%M")
##TODO parse the yaml to get the project_top_level

mkdir -p ${FOLDER}


snakemake  \
--use-conda \
--use-singularity \
--singularity-args "-B /SAN/vyplab:/SAN/vyplab" \
--rerun-triggers mtime \
--jobscript cluster_qsub.sh \
--cluster-config cluster.yaml \
--conda-frontend conda \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--resources gpu=2 \
--nolock \
--rerun-incomplete \
--restart-times 4 \
--latency-wait 100
