#!/usr/bin/env bash

complex=$1

ASTEX_DIR=/ddnB/work/jaydy/dat/astex
WORKING_DIR=/work/jaydy/working/astex_weak/$complex

bin=trace

complex_dir=$ASTEX_DIR/$complex
parameter_dir=/home/jaydy/Workspace/Bitbucket/geauxdock/data/parameters


pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
sdf_file=${complex_dir}/${complex}.sdf
ff_file=${complex_dir}/${complex}-0.8.ff


opt_file=${parameter_dir}/08ff_opt
nor_a_file=${parameter_dir}/08_nor_a
nor_b_file=${parameter_dir}/08_nor_b
para_file=${parameter_dir}/paras

trace_file=$WORKING_DIR/${complex}_grid.csv

mkdir -p $WORKING_DIR
conf_file=$WORKING_DIR/${complex}_new

cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-opt ${opt_file} \
-na ${nor_a_file} \
-nb ${nor_b_file} \
-para ${para_file} \
\
-ns 3000 \
-nc 10 \
-floor_temp 0.044f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
-tr ${trace_file} \
-lc ${conf_file}\
"


echo ${cmd}
${cmd}
