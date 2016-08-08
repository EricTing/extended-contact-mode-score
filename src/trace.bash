#!/bin/bash

#complex=$var1
#echo $complex



gen_conf() {
    complex=$1
    dat=/work/jaydy/dat/astex

    bin=/home/jaydy/Workspace/GitHub/geaxdock_gpu/src/trace
    input_dir=/home/jaydy/Workspace/GitHub/geaxdock_gpu/data
    parameters_dir=${input_dir}/parameters
    complex_dir=${dat}/${complex}

    pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
    sdf_file=${complex_dir}/${complex}.sdf
    ff_file=${complex_dir}/${complex}-0.8.ff

    opt_file=${parameters_dir}/08ff_opt
    nor_a_file=${parameters_dir}/08_nor_a
    nor_b_file=${parameters_dir}/08_nor_b
    para_file=${parameters_dir}/paras

    trace_file=/home/jaydy/Workspace/GitHub/geaxdock_gpu/data/trajectory.sample.out

    output_dir=/work/jaydy/working/astex/rnd_confs/${complex}
    mkdir -p ${output_dir}

    conf_file=${output_dir}/${complex}_new

    # generate new vectors
    python /home/jaydy/Workspace/GitHub/extended-contact-mode-score/src/gen_vec.py > ${trace_file}

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
-r 0.0872f \
-tr ${trace_file} \
-lc ${conf_file}\
"


    echo ${cmd}
    ${cmd}
}

# complex=1a07C1
# gen_conf $complex

for complex in `cat /work/jaydy/dat/astex/total_lst`
do
    echo $complex
    gen_conf $complex
done
