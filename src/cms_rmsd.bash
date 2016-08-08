#!/usr/bin/env bash

rmsd_out=/work/jaydy/working/astex/rmsd.txt
cms_out=/work/jaydy/working/astex/cms.txt


run_rmsd() {
    complex=$1
    dat=/work/jaydy/dat/astex

    bin=/home/jaydy/bin/cms
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

    echo complex $complex >> $rmsd_out

    for num in 0 1 2 3 4 5;
    do
        myconf=${conf_file}_${num}.sdf

        cmd="${bin} -r --lig1 ${myconf} --lig2 ${sdf_file} --prt1 ${pdb_file} --prt2 ${pdb_file}"
        echo $num >> $rmsd_out
        ${cmd} >> $rmsd_out
    done
}


run_cms() {
    complex=$1
    dat=/work/jaydy/dat/astex

    bin=/home/jaydy/bin/cms
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

    echo complex $complex >> $cms_out

    for num in 0 1 2 3 4 5;
    do
        myconf=${conf_file}_${num}.sdf

        cmd="${bin} -c --lig1 ${myconf} --lig2 ${sdf_file} --prt1 ${pdb_file} --prt2 ${pdb_file}"
        echo $num >> $cms_out
        ${cmd} >> $cms_out
    done
}

if [ -f $rmsd_out ] ; then
    rm $rmsd_out
fi

launch_rmsd() {
    for complex in `cat /work/jaydy/dat/astex/total_lst`
    do
        echo $complex
        # run_rmsd $complex
        run_rmsd $complex
    done
}

time launch_rmsd

if [ -f $cms_out ] ; then
    rm $cms_out
fi

launch_cms() {
    for complex in `cat /work/jaydy/dat/astex/total_lst`
    do
        echo $complex
        # run_rmsd $complex
        run_cms $complex
    done
}

time launch_cms
