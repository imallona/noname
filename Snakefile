#!/usr/bin/env snakemake -s
##
## Snakefile to process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op
import os

# include: op.join('src', 'simulate.snmk')
    
configfile: "config.yaml"

include: "src/workflow_functions.py"

## to ease whitelists symlinking
if not op.isabs(config['repo_path']):
    config['repo_path'] = op.join(workflow.basedir, config['repo_path'])

os.makedirs(op.join(config['working_dir'], 'logs'), exist_ok=True)
os.makedirs(op.join(config['working_dir'], 'benchmarks'), exist_ok=True)

print(get_sample_names())



rule all:
    input:
        op.join(config['working_dir'], 'starsolo_wta', 'descriptive_report.html'),
        # 'pbmc_rustody',
        # expand(op.join(config['working_dir'], 'kallisto', '{sample}', 'matrix.ec'),
               # sample = get_sample_names()),
        op.join(config['working_dir'] , 'data', 'mouse_index', 'sampletags', 'SAindex'),
        expand(op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
               sample = get_sample_names())


rule star_index:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        gtf = config['gtf'],
        fa = config['genome']
    output:
        index_path =  op.join(config['working_dir'] , 'data', 'index', 'SAindex')
    threads:
        config['nthreads']
    params:
        simulate = config['simulate'],
        processing_path = op.join(config['working_dir'], 'data'),
        nthreads = config['nthreads'],
        star = config['STAR'],
        sjdbOverhang = config['sjdbOverhang'],
        indexNbases = 4 if config['simulate'] else 14
    log:
        op.join(config['working_dir'], 'logs', 'star_indexing.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'star_indexing.txt')
    shell:
      """
    mkdir -p {params.processing_path}
    # cd {params.processing_path}

    ({params.star} --runThreadN {params.nthreads} \
     --runMode genomeGenerate \
     --sjdbGTFfile {input.gtf} \
     --genomeDir {params.processing_path}/index \
     --genomeSAindexNbases {params.indexNbases} \
     --sjdbOverhang {params.sjdbOverhang} \
     --genomeFastaFiles {input.fa} ) 2> {log}
        """


rule prepare_whitelists:
    # conda:
    #     op.join('envs', 'all_in_one.yaml')
    input:
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample)
    output:
        cb1 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    run:
        sample = wildcards.sample
        symlink_whitelist(sample)


rule starsolo_wta_starsolo:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
        gtf = config['gtf'],
        cb1 = op.join(config['working_dir'], 'starsolo_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    output:
        bam = op.join(config['working_dir'], 'starsolo_wta', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        filtered_barcodes = op.join(config['working_dir'], 'starsolo_wta', '{sample}', 'Solo.out', 'Gene',
                                    'filtered', 'matrix.mtx')
    threads:
        workflow.cores
    log:
        op.join(config['working_dir'], 'logs', 'starsolo_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'star_{sample}.txt')
    params:
        threads = min(10, workflow.cores),
        path = op.join(config['working_dir'], 'starsolo_wta', "{sample}/"),
        index_path = op.join(config['working_dir'] , 'data', 'index'),
        STAR = config['STAR'],
        # num_cells = get_expected_cells_by_name("{sample}"),
        tmp = op.join(config['working_dir'], 'tmp_starsolo_wta_{sample}'),
        maxmem = config['max_mem_mb'] * 1024 * 1024,
        sjdbOverhang = config['sjdbOverhang'],
        soloCellFilter = config['soloCellFilter'],
        soloMultiMappers = config['soloMultiMappers'],
        extraStarSoloArgs = config['extraStarSoloArgs']
    shell:
        """
   rm -rf {params.tmp}
   mkdir -p {params.path} 

   {params.STAR} --runThreadN {params.threads} \
     --genomeDir {params.index_path} \
     --readFilesCommand zcat \
     --outFileNamePrefix {params.path} \
     --readFilesIn  {input.cdna} {input.cbumi}  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist {input.cb1} {input.cb2} {input.cb3} \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter {params.soloCellFilter} \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY \
     --soloCellReadStats Standard \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --soloUMIlen 8 \
     --sjdbGTFfile {input.gtf} \
     --outTmpDir {params.tmp} \
     --sjdbOverhang {params.sjdbOverhang} \
     --limitBAMsortRAM {params.maxmem} \
     --outSAMunmapped Within \
     --soloMultiMappers {params.soloMultiMappers} {params.extraStarSoloArgs}

    rm -rf {params.tmp}
        """

        
# checkpoint retrieve_genome_sizes:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         fa = config['genome']
#     params:
#         faSize = config['faSize']
#     output:
#         op.join(config['working_dir'], 'data', 'chrom.sizes')
#     shell:
#         """
#         {params.faSize} -detailed -tab {input.fa} > {output}
#         """

# ## TODO keep only if generating coverage tracks
# rule index_bam:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         bam = op.join(config['working_dir'], 'starsolo_wta', '{sample}', 'Aligned.sortedByCoord.out.bam')
#     output:
#         bai = op.join(config['working_dir'], 'starsolo_wta', '{sample}',
#                       'Aligned.sortedByCoord.out.bam.bai')
#     threads: workflow.cores
#     shell:
#         """
#         samtools index -@ {threads} {input.bam}     
#         """

## yes the log is considered an output - to pass as a flag
## R_LIBS are conda's if run in conda, but /home/rock/R_LIBs if run in docker, and user's if run directly
rule install_r_deps:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        script = op.join(config['repo_path'], 'src', 'installs.R')
    output:
        log = op.join(config['working_dir'], 'log', 'installs.log')
    params:
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin'],
    log:
        op.join(config['working_dir'], 'logs', 'r_install.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_install.txt')
    shell:
        """
        mkdir -p {params.log_path}

        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} &> {log}
         """
        
rule generate_sce:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        wta_filtered = op.join(config['working_dir'], 'starsolo_wta', '{sample}', 'Solo.out',
                               'Gene', 'filtered', 'matrix.mtx'),
        # gtf = config['gtf'],
        script = op.join(config['repo_path'], 'src', 'generate_sce_object.R'),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'starsolo_wta', '{sample}', '{sample}_sce.rds')
    params:
        align_path = op.join(config['working_dir'], 'starsolo_wta'),
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin']
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}.txt')
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --output_fn {output.sce} &> {log}
        """

                
# rule generate_sce_tasseq:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         wta_filtered = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out',
#                                'Gene', 'filtered', 'matrix.mtx'),
#         # gtf = config['gtf'],
#         script = op.join(config['repo_path'], 'src', 'generate_sce_object.R'),
#         installs = op.join(config['working_dir'], 'log', 'installs.log')
#     output:
#         sce = op.join(config['working_dir'], 'starsolo_wta', '{sample}', '{sample}_tasseq_sce.rds')
#     params:
#         align_path = op.join(config['working_dir'], 'tasseq'),
#         working_dir = config['working_dir'],
#         sample = "{wildcards.sample}",
#         Rbin = config['Rbin']
#     shell:
#         """
#         {params.Rbin} -q --no-save --no-restore --slave \
#              -f {input.script} --args \
#              --sample {wildcards.sample} \
#              --working_dir {params.working_dir} \
#              --output_fn {output.sce}
#         """


rule render_descriptive_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        # mapping_report = op.join(config['working_dir'], 'multimodal', 'mapping_summary.txt'),
        # gtf = config['gtf'],
        script = op.join(config['repo_path'], 'src', 'generate_descriptive_singlecell_report.Rmd'),
        sces = expand(op.join(config['working_dir'], 'starsolo_wta', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names()),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        html = op.join(config['working_dir'], 'starsolo_wta', 'descriptive_report.html')
        # cache = temp(op.join(config['repo_path'], 'process_sce_objects_cache')),
        # cached_files = temp(op.join(config['repo_path'], 'process_sce_objects_files'))
    log:
        op.join(config['working_dir'], 'logs', 'descriptive_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'descriptive_report.log')
    params:
        path = op.join(config['working_dir'], 'starsolo_wta'),
        working_dir = op.join(config['working_dir']),
        sample = "{wildcards.sample}",
        Rbin = config['Rbin']
    shell:
        """
        cd {params.working_dir}
        mkdir -p {params.path}
        {params.Rbin} --vanilla -e 'rmarkdown::render(\"{input.script}\", 
          output_file = \"{output.html}\", 
          params = list(path = \"{params.path}\"))' &> {log}
        """

rule install_rustody:
    conda:
        op.join('envs', 'all_in_one.yaml')
    output:
        op.join('soft', 'Rustody', 'target', '.rustc_info.json')
    log:
        op.join(config['working_dir'], 'logs', 'rustody_install.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'rustody_install.txt')
    shell:
        """
        mkdir -p soft
        cd soft
        curl https://sh.rustup.rs -sSf | sh
        source "$HOME/.cargo/env"
        git clone https://github.com/stela2502/Rustody --depth 1
        cd Rustody
        cargo build --release 2> {log}
        """
        
rule rustody_CAUTIONMAXREADS:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        rustody = op.join('soft', 'Rustody', 'target', '.rustc_info.json'),
        transcriptome = config['transcriptome'],
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cb_umi = lambda wildcards: get_cbumi_by_name(wildcards.sample)        
    output:
        flag = '{sample}_rustody'
    threads:
        workflow.cores
    params:
        whitelist = lambda wildcards: get_barcode_whitelist_by_name(wildcards.sample),
        species = lambda wildcards: get_species_by_name(wildcards.sample),
        rustody_path = op.join('soft', 'Rustody', 'target', 'release')
    log:
        op.join(config['working_dir'], 'logs', 'rustody_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'rustody_{sample}.txt')
    shell:
        """
        source "$HOME/.cargo/env"
        if [ {params.whitelist} == '384x3' ]; then
           wl='v2.384'
        elif [ {params}.whitelist == '96x3' ]; then
           wl='v2.96'
        else
           wl='error_unknown_whitelist_spec'
        fi

        export PATH="{params.rustody_path}:"$PATH
        quantify_rhapsody_multi \
           --version "$wl" \
           --specie {params.species} \
           --reads {input.cb_umi} \
           --outpath ~/Rustody_pdgfra \
           --num-threads {threads} \
           --file {input.cdna} \
           --expression {input.transcriptome} \
           --min-umi 1 \
           --max-reads 100000 $ > {log}

        touch {output.flag}
        """

rule standardize_cb_umis_cutadapt:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        cb_umi = lambda wildcards: get_cbumi_by_name(wildcards.sample)
    output:
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz")
    params:
        path = op.join(config['working_dir'], 'data', 'fastq')
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.path}
        pigz -p {threads} --decompress --stdout {input.cb_umi} | \
           cutadapt -g "NNNNNNNNNGTGANNNNNNNNNGACANNNNNNNNNNNNNNNNN;min_overlap=43;noindels" \
           --action=crop \
           --discard-untrimmed \
           --cores {threads} \
           -e 0 \
           - | cut -c1-9,14-22,27- | pigz -p {threads} > {output.standardized_cb_umi}
        """

## conda recipe is broken        
rule kallisto_index:
    # conda:
    #     op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome']
    params:
        index_name = op.join('data', 'index', 'kallisto_index')
    output:
        op.join(config['working_dir'] , 'data', 'index', 'kallisto')
    threads:
        workflow.cores
    log:
        op.join(config['working_dir'], 'logs', 'kallisto_index.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'kallisto_index.txt')
    shell:
        """
        kallisto index --threads={threads} \
           -i={params.index_name} {input.transcriptome} &> {log}
        """
## conda recipe is broken 
rule kallisto_bus:
    # conda:
    #     op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome'],
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
        kallisto_index = op.join(config['working_dir'] , 'data', 'index', 'kallisto')
    output:
        op.join(config['working_dir'], 'kallisto', '{sample}', 'matrix.ec')
    params:
        output_dir = op.join(config['working_dir'], 'kallisto', '{sample}')
    threads: workflow.cores
    log:
        op.join(config['working_dir'], 'logs', '{sample}_kallisto_bus.log')
    log:
        op.join(config['working_dir'], 'benchmarks', '{sample}_kallisto_bus.txt')
    shell:
        """
        kallisto bus --index {input.index} \
            --output-dir {params.output_dir} \
            -x '0,0,27:0,27,35:1,0,0' \
            -t {threads} \
            {input.standardized_cb_umi} {input.cdna} &> {log}
        """

## from https://github.com/imallona/rock_roi_paper/blob/imallona/03_leukemia/02_sampletags_again.sh

for sample in get_sample_names():
    species = get_species_by_name(name = sample)
    rule star_index_sampletags:
        conda:
            op.join('envs', 'all_in_one.yaml')
        input:
            fa = op.join('data', 'sampletags', species + '_sampletags.fa')
        output:
            op.join(config['working_dir'] , 'data', species + '_index', 'sampletags', 'SAindex')
        threads:
            workflow.cores
        params:
            output_dir = op.join(config['working_dir'], 'data', species + '_index', 'sampletags')
        log:
            op.join(config['working_dir'], 'logs', species + '_sampletags_index.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', species + '_sampletags__index.txt')
        shell:
            """
            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeSAindexNbases 2 \
            --genomeDir {params.output_dir} \
            --genomeFastaFiles {input.fa} &> {log}
            """

# rule extract_unmapped_startsolo_wta:
    
        
## rule align_starsolo_sampletags:
    # input:
    # output:
    # params:
    # log:
    # benchmark:
    # shell:
    
# # https://github.com/s-shichino1989/TASSeq_EnhancedBeads/blob/e48fd2c2fd5a23d622f03e206b8fbe87772fd57f/shell_scripts/Rhapsody_STARsolo.sh#L18
# rule starsolo_wta_tasseq_style:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
#         cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
#         index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
#         gtf = config['gtf'],
#         cb1 = op.join(config['working_dir'], 'starsolo_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
#         cb2 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
#         cb3 = op.join(config['working_dir'], 'starsolo_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
#     output:
#         bam = op.join(config['working_dir'], 'tasseq', '{sample}', 'Aligned.sortedByCoord.out.bam'),
#         filtered_barcodes = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out', 'Gene',
#                                     'filtered', 'barcodes.tsv'),
#         filtered_counts = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out', 'Gene',
#                                   'filtered', 'matrix.mtx')
#     threads: workflow.cores
#     params:
#         threads = min(10, workflow.cores),
#         path = op.join(config['working_dir'], 'tasseq', "{sample}/"),
#         index_path = op.join(config['working_dir'] , 'data', 'index'),
#         STAR = config['STAR'],
#         # num_cells = get_expected_cells_by_name("{sample}"),
#         tmp = op.join(config['working_dir'], 'tmp_tasseq_{sample}'),
#         maxmem = config['max_mem_mb'] * 1024 * 1024,
#         sjdbOverhang = config['sjdbOverhang'],
#         soloCellFilter = config['soloCellFilter']
#     shell:
#         """
#         rm -rf {params.tmp}
#         mkdir -p {params.path} 

#         {params.STAR} --runThreadN {params.threads} \
#           --genomeDir {params.index_path} \
#         --readFilesIn {input.cdna} {input.cbumi} \
#         --outFileNamePrefix {params.path} \
#         --readFilesCommand zcat \
#         --clipAdapterType CellRanger4 \
#         --outSAMtype BAM SortedByCoordinate \
#         --outBAMsortingThreadN {params.threads} \
#         --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
#         --outSAMunmapped Within \
#         --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
#         --outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
#         --soloCellFilter {params.soloCellFilter} \
#         --soloUMIdedup Exact \
#         --soloMultiMappers Rescue \
#         --soloFeatures Gene GeneFull \
#         --soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
#         --soloCBmatchWLtype EditDist_2 \
#         --soloCBwhitelist {input.cb1} {input.cb2} {input.cb3} \
#         --soloType CB_UMI_Complex \
#         --soloUMIlen 8 \
#         --soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
#         --soloUMIposition 3_10_3_17 
#     rm -rf {params.tmp}
#         """
