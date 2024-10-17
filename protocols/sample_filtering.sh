set -e
set -u

ml ${plink2Version}
ml ${RPlusVersion}
export R_LIBS_USER=${rLibsPath}

mkdir -p $(dirname ${sampleListPrefix})

plink2 --bfile ${concatenatedGenotypesOutputDir}/plink_dataset_to_qc --out 'plink2_sample_qc' --missing sample-only --het

Rscript ${pipelineRoot}/scripts/qc_autosomes.R \
--sample-missingness plink2_sample_qc.smiss \
--heterozygosity plink2_sample_qc.het \
--out-prefix qc_out

plink2 --bfile ${concatenatedGenotypesOutputDir}/plink_dataset_to_qc --keep qc_out.samples_passed_qc.txt \
--make-bed --out ${concatenatedGenotypesOutputDir}/chr_all

awk 'BEGIN{FS="\t"; OFS=FS}{print $2}' ${concatenatedGenotypesOutputDir}/chr_all.fam > ${sampleListPrefix}.samples.txt
