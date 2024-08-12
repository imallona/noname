#!/usr/bin/env snakemake -s
##
## Snakefile to process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op

# include: op.join('src', 'simulate.snmk')
    
configfile: "config.yaml"

include: "src/workflow_functions.py"

## to ease whitelists symlinking
if not op.isabs(config['repo_path']):
    config['repo_path'] = op.join(workflow.basedir, config['repo_path'])

print(get_sample_names())

rule all:
    input:
        # expand(op.join(config['working_dir'], 'align_wta', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        #                sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'align_wta', '{sample}', '{sample}_sce.rds'),
        #        sample = get_sample_names()),
        op.join(config['working_dir'], 'align_wta', 'descriptive_report.html'),
        'pbmc_flag'
        

rule index:
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
        op.join(config['working_dir'], 'data', 'indexing.log')
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
        # r1 = op.join(config['working_dir'], 'data', "{sample}", 'r1.fq.gz'),
        # r2 = op.join(config['working_dir'], 'data', "{sample}", 'r2.fq.gz')
    output:
        cb1 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    run:
        sample = wildcards.sample
        symlink_whitelist(sample)


rule align_wta:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
        gtf = config['gtf'],
        cb1 = op.join(config['working_dir'], 'align_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    output:
        bam = op.join(config['working_dir'], 'align_wta', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                    'filtered', 'barcodes.tsv')
    threads: workflow.cores
    params:
        threads = min(10, workflow.cores),
        path = op.join(config['working_dir'], 'align_wta', "{sample}/"),
        index_path = op.join(config['working_dir'] , 'data', 'index'),
        STAR = config['STAR'],
        # num_cells = get_expected_cells_by_name("{sample}"),
        tmp = op.join(config['working_dir'], 'tmp_align_wta_{sample}'),
        maxmem = config['max_mem_mb'] * 1024 * 1024,
        sjdbOverhang = config['sjdbOverhang'],
        soloCellFilter = config['soloCellFilter']
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
     --limitBAMsortRAM {params.maxmem}

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

## TODO keep only if generating coverage tracks
rule index_bam:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        bam = op.join(config['working_dir'], 'align_wta', '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        bai = op.join(config['working_dir'], 'align_wta', '{sample}',
                      'Aligned.sortedByCoord.out.bam.bai')
    threads: workflow.cores
    shell:
        """
        samtools index -@ {threads} {input.bam}     
        """

# rule generate_tracks:
    
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
        log_path = op.join(config['working_dir'], 'log')
    shell:
        """
        mkdir -p {params.log_path}

        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} &> {output.log}
#         """
        
rule generate_sce:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        wta_filtered = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out',
                               'Gene', 'filtered', 'matrix.mtx'),
        # gtf = config['gtf'],
        script = op.join(config['repo_path'], 'src', 'generate_sce_object.R'),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'align_wta', '{sample}', '{sample}_sce.rds')
    params:
        align_path = op.join(config['working_dir'], 'align_wta'),
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin']
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --output_fn {output.sce}
        """


rule render_descriptive_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        # mapping_report = op.join(config['working_dir'], 'multimodal', 'mapping_summary.txt'),
        # gtf = config['gtf'],
        script = op.join(config['repo_path'], 'src', 'generate_descriptive_singlecell_report.Rmd'),
        sces = expand(op.join(config['working_dir'], 'align_wta', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names()),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        html = op.join(config['working_dir'], 'align_wta', 'descriptive_report.html')
        # cache = temp(op.join(config['repo_path'], 'process_sce_objects_cache')),
        # cached_files = temp(op.join(config['repo_path'], 'process_sce_objects_files'))
    log: op.join(config['working_dir'], 'align_wta', 'descriptive_report.log')
    params:
        path = op.join(config['working_dir'], 'align_wta'),
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
        'rustody_install.log'
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
        
rule rustody:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        rustody = op.join('soft', 'Rustody', 'target', '.rustc_info.json'),
        transcriptome = config['transcriptome'],
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cb_umi = lambda wildcards: get_cbumi_by_name(wildcards.sample)        
    output:
        flag = '{sample}_flag'
    threads: workflow.cores
    params:
        whitelist = lambda wildcards: get_barcode_whitelist_by_name(wildcards.sample),
        species = lambda wildcards: get_species_by_name(wildcards.sample),
        rustody_path = op.join('soft', 'Rustody', 'target', 'release')
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
           --min-umi 100

        touch {output.flag}
        """

