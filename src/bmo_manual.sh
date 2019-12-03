#!/bin/bash

# Config
bmo_dir="."
sample_name="example_buenrostro_rep1"
bam="examples/${sample_name}.subsampled.bam"
peaks="examples/example_peaks.bed.gz"
motif_list="examples/motifs.txt"
motif_dir="examples"
motif_ext=".bed.gz"
cores="4"


# Run
raw_signals_dir="work/raw_signals/${sample_name}/all_regions"
outside_peaks_dir="work/raw_signals/${sample_name}/outside_peaks"
cooccur_dir="work/raw_signals/${sample_name}/co_occurring"
nb_dir="work/negative_binomials/${sample_name}"
bmo_outdir="work/bmo/${sample_name}"

mkdir -p ${raw_signals_dir} ${outside_peaks_dir} ${cooccur_dir} \
    ${nb_dir} ${bmo_outdir} 


# Pipeline

echo "# measure_raw_signal"
echo "# drmr:job processors=4 processor_memory=1gb time_limit=15:00"
while read motif
do
    motif_file="${motif_dir}/${motif}${motif_ext}"
    raw_signals="${raw_signals_dir}/${motif}.bed"
    
    echo "${bmo_dir}/src/measureRawSignal.py -p ${cores} -b ${bam} \
        -m ${motif_file} > ${raw_signals}"
done < ${motif_list}


echo -e "\n\n# motifs_outside_peaks"
echo "# drmr:job processors=1 processor_memory=1gb time_limit=15:00"
while read motif
do
    raw_signals="${raw_signals_dir}/${motif}.bed"
    outside_peaks="${outside_peaks_dir}/${motif}.bed"
    
    echo "intersectBed -v -a ${raw_signals} -b ${peaks} > ${outside_peaks}"
done < ${motif_list}


echo -e "\n\n# count_overlapping_motifs"
echo "# drmr:job processors=1 processor_memory=1gb time_limit=15:00"
while read motif
do
    raw_signals="${raw_signals_dir}/${motif}.bed"
    cooccur="${cooccur_dir}/${motif}.bed"

    echo "${bmo_dir}/src/count_co-occuring_motifs.sh ${raw_signals} 100 \
        ${bmo_dir} > ${cooccur}"
done < ${motif_list}


echo -e "\n\n# fit_nbinoms"
echo "# drmr:job processors=1 processor_memory=5gb time_limit=30:00"
while read motif
do
    nb_handle="${nb_dir}/${motif}"
    nb_out="${nb_handle}.bed.gz"
    echo "Rscript ${bmo_dir}/src/rawSigNBModel.R -f ${motif}.bed \
            --dir1 ${raw_signals_dir}/ --dir2 ${outside_peaks_dir}/ \
            -o ${nb_handle} -c 7 --writeBed"
done < ${motif_list}


echo -e "\n\n# BMO"
echo "# drmr:job processors=1 processor_memory=5gb time_limit=15:00"
while read motif
do
    cooccur="${cooccur_dir}/${motif}.bed"
    nb_handle="${nb_dir}/${motif}"
    nb_out="${nb_handle}.bed.gz"

    echo "Rscript ${bmo_dir}/src/bmo.R --f1 ${nb_out} --f2 ${cooccur} \
        -o ${bmo_outdir} -n ${motif} -p 1"
done < ${motif_list}

