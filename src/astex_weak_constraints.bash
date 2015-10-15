#!/usr/bin/env bash


run() {
    local complex=$1
    ASTEX_DIR=/ddnB/work/jaydy/dat/astex
    WORKING_DIR=/work/jaydy/working/astex_weak/$complex

    bin=~/bin/dock

    complex_dir=$ASTEX_DIR/$complex
    parameter_dir=/home/jaydy/Workspace/Bitbucket/geauxdock/data/parameters


    pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
    sdf_file=${complex_dir}/${complex}.sdf
    ff_file=${complex_dir}/${complex}-0.8.ff


    para_file=${parameter_dir}/paras


    mkdir -p $WORKING_DIR
    csv_file=$WORKING_DIR/${complex}.csv

    cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-para ${para_file} \
\
-nc 10 \
-floor_temp 0.004f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
-csv $csv_file \
"


    echo ${cmd}
    ${cmd}
}

for c in `cat /home/jaydy/Workspace/Bitbucket/xcms/dat/astex.lst`
do
    run $c
done
