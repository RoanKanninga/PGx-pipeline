nextflow.enable.dsl=2

def helpMessage() {
    log.info """
    eqtlgenimpute v${workflow.manifest.version}
    Usage: nextflow run main.nf --bfile PREFIX --outdir OUTDIR --output_name NAME
    """.stripIndent()
}
if (params.help) { helpMessage(); exit 0 }

process crossmap {
    tag 'crossmap'
    input:
    tuple path(bed), path(bim), path(fam)
    path chain_file
    output:
    tuple path('crossmapped_plink.bed'), path('crossmapped_plink.bim'), path('crossmapped_plink.fam')
    script:
    """
    awk '{print \\$1,\\$4,\\$4+1,\\$2,\\$5,\\$6,\\$2 "___" \\$5 "___" \\$6}' ${bim} > crossmap_input.bed
    CrossMap.py bed ${chain_file} crossmap_input.bed crossmap_output.bed
    awk '{print \\$7}' crossmap_input.bed | sort > input_ids.txt
    awk '{print \\$7}' crossmap_output.bed | sort > output_ids.txt
    comm -23 input_ids.txt output_ids.txt | awk '{split(\\$0,a,"___"); print a[1]}' > excluded_ids.txt
    plink2 --bfile ${bed.baseName} --exclude excluded_ids.txt --make-bed --output-chr MT --out crossmapped_plink --keep-allele-order --update-map crossmap_output.bed 2 4
    """
}

process sort_bed {
    input: tuple path(bed), path(bim), path(fam)
    output: tuple path('sorted.bed'), path('sorted.bim'), path('sorted.fam')
    script:
    """
    plink2 --bfile ${bed.baseName} --make-bed --output-chr MT --out sorted
    """
}

process plink_to_vcf {
    input: tuple path(bed), path(bim), path(fam)
    output: path('sorted_hg38.vcf')
    script:
    """
    plink2 --bfile ${bed.baseName} --export vcf-4.2 id-paste=iid --chr 1-22 --out sorted_hg38
    """
}

process fix_indels {
    input:
    path vcf
    path fasta
    output: path('sorted_genotypes_hg38_indels.vcf')
    script:
    """
    samtools faidx ${fasta}
    fixvcf.py -f ${fasta} -o sorted_genotypes_hg38_indels.vcf ${vcf}
    """
}

process vcf_fixref_hg38 {
    input:
    path input_vcf
    path fasta
    tuple path(ref_vcf), path(ref_index)
    output: path('fixref_hg38.vcf.gz')
    script:
    """
    bcftools sort ${input_vcf} -Oz -o input.vcf.gz
    bcftools index input.vcf.gz
    bcftools +fixref input.vcf.gz -- -f ${fasta} -i ${ref_vcf} |
      bcftools norm --check-ref x -f ${fasta} |
      bcftools sort -Oz -o fixref_hg38.vcf.gz
    """
}

process filter_preimpute_vcf {
    publishDir "${params.outdir}/preimpute/", mode: 'copy'
    input: path input_vcf
    output:
    tuple path('filtered.vcf.gz'), path('filtered.vcf.gz.csi')
    script:
    """
    bcftools sort ${input_vcf} -Oz -o input.vcf.gz
    bcftools index input.vcf.gz
    bcftools +fill-tags input.vcf.gz -Oz -o tagged.vcf.gz
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05' tagged.vcf.gz |
      bcftools norm -d all |
      bcftools norm -m+any |
      bcftools view -m2 -M2 -Oz -o filtered.vcf.gz
    bcftools index --csi filtered.vcf.gz
    """
}

process calculate_missingness {
    publishDir "${params.outdir}/preimpute/", mode: 'copy'
    input: tuple path(vcf), path(index)
    output: path('genotypes.imiss')
    script:
    """
    vcftools --gzvcf ${vcf} --missing-indv --out genotypes
    """
}

process split_by_chr {
    publishDir "${params.outdir}/preimpute/split_chr", mode: 'copy'
    input:
    tuple path(vcf), path(index)
    val chr
    output: tuple val(chr), path("chr_${chr}.vcf.gz")
    script:
    """
    bcftools view -r ${chr} ${vcf} -Oz -o chr_${chr}.vcf.gz
    """
}

process eagle_prephasing {
    publishDir "${params.outdir}/phasing", mode: 'copy'
    input:
    tuple val(chromosome), path(vcf)
    path genetic_map
    path phasing_reference
    output:
    tuple val(chromosome), path("chr${chromosome}.phased.vcf.gz"), path("chr${chromosome}.phased.vcf.gz.tbi")
    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} --vcfRef=chr${chromosome}.bcf \\
      --geneticMapFile=${genetic_map} --chrom=${chromosome} \\
      --outPrefix=chr${chromosome}.phased --numThreads=${task.cpus} --noImpMissing
    tabix chr${chromosome}.phased.vcf.gz
    """
}

process minimac_imputation {
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: '*.dose.vcf.gz'
    input:
    tuple val(chromosome), path(vcf), path(index), val(chrom), val(start), val(end), val(name)
    path imputation_reference
    val flank_size
    output:
    tuple val(chromosome), val(start), val(end), val(name), path("range_${chromosome}_${start}-${end}_${name}.dose.vcf.gz")
    script:
    """
    minimac4 --refHaps chr${chromosome}.m3vcf.gz --haps ${vcf} --chr ${chromosome} \\
      --start ${start} --end ${end} --window ${flank_size} --prefix imputed \\
      --format GT,GP,HDS --noPhoneHome
    tabix imputed.dose.vcf.gz
    bcftools concat ${vcf} imputed.dose.vcf.gz --regions ${chromosome}:${start}-${end} \\
      --remove-duplicates --allow-overlaps \\
      -Oz -o range_${chromosome}_${start}-${end}_${name}.dose.vcf.gz
    """
}

process index_imputation {
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: '*.tbi'
    input: tuple val(chromosome), val(start), val(end), val(name), path(vcf)
    output: tuple val(chromosome), val(start), val(end), val(name), path(vcf), path('*.tbi')
    script: "tabix ${vcf}"
}

process annotate_imputation {
    publishDir "${params.outdir}/annotated/", mode: 'copy', pattern: '*.annotated.bar.vcf.gz*'
    input:
    tuple val(chromosome), val(start), val(end), val(name), path(vcf), path(index)
    path annotation_vcf
    output:
    tuple val(chromosome), val(start), val(end), val(name),
          path("range_${chromosome}_${start}-${end}_${name}.annotated.bar.vcf.gz"),
          path("range_${chromosome}_${start}-${end}_${name}.annotated.bar.vcf.gz.tbi")
    script:
    """
    tabix ${annotation_vcf}
    bcftools annotate ${vcf} -x 'INFO/AF' --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o removed_af_info.vcf.gz
    tabix removed_af_info.vcf.gz
    bcftools annotate --annotations ${annotation_vcf} --columns 'INFO,ID' -Ov removed_af_info.vcf.gz > removed_af_info.TMP.vcf
    bcftools +setGT removed_af_info.TMP.vcf -- -t . -n p | bgzip -c > range_${chromosome}_${start}-${end}_${name}.annotated.bar.vcf.gz
    tabix range_${chromosome}_${start}-${end}_${name}.annotated.bar.vcf.gz
    """
}

process split_target_dataset {
    publishDir "${params.outdir}/target/", mode: 'copy', pattern: '*.target.vcf.gz*'
    input:
    tuple val(chromosome), val(start), val(end), val(name), path(vcf), path(index)
    path target_samples
    output:
    tuple val(chromosome), val(start), val(end), val(name),
          path("range_${chromosome}_${start}-${end}_${name}.annotated.bar.target.vcf.gz"),
          path("range_${chromosome}_${start}-${end}_${name}.annotated.bar.target.vcf.gz.tbi")
    script:
    """
    bcftools view --samples-file ${target_samples} --force-samples ${vcf} -Ov |
      bcftools +fill-tags -Ov -- -t AF |
      bgzip -c > range_${chromosome}_${start}-${end}_${name}.annotated.bar.target.vcf.gz
    tabix range_${chromosome}_${start}-${end}_${name}.annotated.bar.target.vcf.gz
    """
}

workflow {
    // Input channels
    bfile_ch = Channel.value(tuple(file("${params.bfile}.bed"), file("${params.bfile}.bim"), file("${params.bfile}.fam")))
    chain_ch = Channel.fromPath(params.chain_file).ifEmpty { error "Chain file not found: ${params.chain_file}" }.collect()
    fasta_ch = Channel.fromPath(params.target_ref).ifEmpty { error "Target reference not found: ${params.target_ref}" }.collect()
    ref_panel_ch = Channel.fromPath("${params.ref_panel_hg38}*").ifEmpty { error "Reference panel not found: ${params.ref_panel_hg38}" }.collect()
    map_ch = Channel.fromPath(params.eagle_genetic_map).ifEmpty { error "Eagle map not found: ${params.eagle_genetic_map}" }.collect()
    eagle_ref_ch = Channel.fromPath("${params.eagle_phasing_reference}*").ifEmpty { error "Eagle reference not found: ${params.eagle_phasing_reference}" }.collect()
    minimac_ref_ch = Channel.fromPath("${params.minimac_imputation_reference}*").ifEmpty { error "Minimac reference not found: ${params.minimac_imputation_reference}" }.collect()
    annotation_ch = Channel.fromPath(params.annotation_vcf_file).collect()
    ranges_ch = Channel.fromPath(params.range_bed_hg38)
        .splitCsv(header: ['chrom','start','end','name'], sep: '\\t')
        .map { row -> tuple(row.chrom, row.start, row.end, row.name) }
    chr_ch = ranges_ch.map { it[0] }.unique()

    crossmap(bfile_ch, chain_ch)
    sort_bed(crossmap.out)
    plink_to_vcf(sort_bed.out)
    fix_indels(plink_to_vcf.out, fasta_ch)
    vcf_fixref_hg38(fix_indels.out, fasta_ch, ref_panel_ch)
    filter_preimpute_vcf(vcf_fixref_hg38.out)
    calculate_missingness(filter_preimpute_vcf.out)
    split_by_chr(filter_preimpute_vcf.out, chr_ch)
    eagle_prephasing(split_by_chr.out, map_ch, eagle_ref_ch)

    // Preserve chromosome/range association explicitly instead of DSL1 cross().
    phased_with_ranges = eagle_prephasing.out
        .combine(ranges_ch, by: 0)
        .map { phased, range -> tuple(phased[0], phased[1], phased[2], range[0], range[1], range[2], range[3]) }
    minimac_imputation(phased_with_ranges, minimac_ref_ch, Channel.value(params.imputation_flank_size))
    index_imputation(minimac_imputation.out)
    annotate_imputation(index_imputation.out, annotation_ch)
    target_samples_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .filter { it.analysis == 'diagnostics' }
        .map { it.Sample_ID }
        .collectFile(name: 'target_samples.txt', newLine: true)
    split_target_dataset(annotate_imputation.out, target_samples_ch)
}
