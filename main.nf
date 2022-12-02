process REFERENCE_INDEX {
    
    //publishDir "$params.results/reference_index" 
    input:
        path ref
    output:
        tuple path(ref), path("*.fai")
    shell:
        '''
        samtools faidx !{ref}
        '''

}

process BAM_INDEX {
    
    //publishDir "$params.results/bam_index" 
    input:
        path bam
    output:
        tuple path(bam), path("*.bai")
    shell:
        '''
        samtools index -@!{task.cpus} !{bam}
        '''

}

process BREAKDANCER {
    
    publishDir "$params.results/breakdancer" 
    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        val map
        val filt
    output:
        file "*breakdancer.vcf"
        
    shell:
        def my_map = map ? map : 20
        def my_filt = filt ? filt : 20
        prefix = "${bam.baseName}.breakdancer"
        '''
        bam2cfg.pl -q !{my_map} -g !{bam} > "!{prefix}.cfg"
        breakdancer-max -q !{my_map} -y !{my_filt} -h "!{prefix}.cfg" > "!{prefix}.ctx"
        breakdancertovcf.py -o "!{prefix}.vcf" !{ref} "!{prefix}.ctx"
        
        '''

}

process DELLY {
    
    publishDir "$params.results" 
    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed 
    output:
        file "*delly.vcf"
    shell:
        def my_bed = bed ? " -x $bed" : "" 
        outfile = "${bam.simpleName}.delly.vcf"
        '''
        delly call!{my_bed} -t ALL -g !{ref} -o delly_output.bcf !{bam}
        bcftools view -Ov -o !{outfile} delly_output.bcf 
        '''

}

process LUMPY {
    
    publishDir "$params.results" 
    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed
    output:
        file "*lumpy.vcf"
    shell:
        def my_bed = bed ? " --exclude $bed" : "" 
        prefix = bam.simpleName
        outfile = "${bam.simpleName}.lumpy.vcf.gz"
        '''
        smoove call -x --genotype --name !{prefix} --outdir . -f !{ref} --processes !{task.cpus}!{my_bed} !{bam}
        mv *genotyped.vcf.gz !{outfile}
        bgzip -d !{outfile} 
        
        '''

}

process MANTA {
    
    publishDir "$params.results"
    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed 
    output:
        file "*manta.vcf"
    shell:
        def my_bed = bed ? " --callRegions $bed" : ""  // To update with include regions not finish
        outfile = "${bam.simpleName}.manta.vcf"
        '''
        MEMORY="!{task.memory}"
        MEM=${MEMORY/ GB/}
        configManta.py --bam !{bam} --referenceFasta !{ref} --runDir ./
        ./runWorkflow.py -j !{task.cpus} -g $MEM
        cp results/variants/diploidSV.vcf.gz .
        bgzip -d diploidSV.vcf.gz
        mv diploidSV.vcf !{outfile}
        '''

}

process PINDEL_SINGLE {
    
    publishDir "$params.results/pindel" 
    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        each chr
    output:
        file "*.vcf"
    shell:
    fileName = bam.name
    fileSimpleName = bam.simpleName
    outfile = "${fileSimpleName}.pindel.${chr}.vcf"
        '''
        echo -e !{fileName}'\t'350'\t'!{fileSimpleName} > pindel_config.txt
        pindel -T !{task.cpus} -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}
        pindel2vcf -P !{fileSimpleName} -r !{ref} -R GRCH37 -d `date +'%m/%d/%Y'` -v !{outfile}

        '''

}

process MERGE_PINDEL_SINGLE {
    
    publishDir "$params.results" 
    input:
        path vcfs
    output:
        file "*.pindel.vcf"
    shell:
        '''
        for f in $(ls *.vcf); do bgzip $f;tabix -p vcf "${f}.gz";done
        bcftools concat -Ov -o "${f%%.*}.pindel.unsorted.vcf" *.vcf.gz
        bcftools sort -Ov -o "${f%%.*}.pindel.vcf" "${f%%.*}.pindel.unsorted.vcf"

        '''

}

process TARDIS_PREP {
    
    publishDir "$params.results/tardis_prep" 
    input:
        tuple path(bam),path(bam_index)

    output:
        tuple path("*.markdup.bam"), path("*.markdup.bam.bai")
    shell:
        outfile = "${bam.baseName}.markdup.bam"
        '''
        mkdir tmp_sam
        sambamba markdup -r -t !{task.cpus} --tmpdir=tmp_sam !{bam} !{outfile}
        samtools index -@!{task.cpus} !{outfile}

        '''

}

process TARDIS {
    
    publishDir "$params.results" 
    input:
        tuple path(bam), path(bam_index)
        tuple path(ref),path(ref_index)
        path sonic
        path bed 
    output:
        file "*.tardis.vcf"
    shell:
        def my_bed = bed ? " --gaps $bed" : "" 
        outfile = "${bam.baseName}"
        '''
        tardis -i !{bam} --ref !{ref} --sonic !{sonic} --out !{outfile}!{my_bed}
        cat "!{outfile}.vcf" | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.tardis.vcf"

        '''

}

workflow {
    chromosomes = Channel.of(1..22,"X","Y")

    REFERENCE_INDEX (params.reference)
    BAM_INDEX (params.input)
    BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out)
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, params.exclude_bed)
    MANTA (BAM_INDEX.out, REFERENCE_INDEX.out)
    PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, chromosomes)
    TARDIS_PREP (BAM_INDEX.out)
    TARDIS(TARDIS_PREP.out, REFERENCE_INDEX.out, params.sonic_file)
    MERGE_PINDEL_SINGLE(PINDEL_SINGLE.out.collect())
}
