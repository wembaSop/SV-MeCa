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

process INCLUDE_REGIONS {

    input:
        tuple path(bam), path(bai) 
        path bed

    output:
        path "include.bed"

    shell:

        '''
        samtools view -H !{bam}| grep "^@SQ"| awk -F '[\t:]' '{print $3"\t"0"\t"$5}'> full.bed
        bedtools subtract -a full.bed -b !{bed} > include.bed

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

        prefix = "${bam.baseName}.breakdancer"

        '''
        bam2cfg.pl -q !{map} -g !{bam} > "!{prefix}.cfg"
        breakdancer-max -q !{map} -y !{filt} -h "!{prefix}.cfg" > "!{prefix}.ctx"
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
        path bed //regions to include

    output:
        file "*manta.vcf"

    shell:

        def my_bed = bed ? " --callRegions $bed" : ""  
        outfile = "${bam.simpleName}.manta.vcf"

        '''
        MEMORY="!{task.memory}"
        MEM=${MEMORY/ GB/}
        configManta.py --bam !{bam} --referenceFasta !{ref} --runDir ./!{my_bed}
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
        path bed
        val insert
        each chr

    output:
        file "*.vcf"

    shell:

        def my_bed = bed ? " --exclude $bed" : ""
        fileName = bam.name
        fileSimpleName = bam.simpleName
        outfile = "${fileSimpleName}.pindel.${chr}.vcf"

        '''
        echo -e !{fileName}'\t'!{insert}'\t'!{fileSimpleName} > pindel_config.txt
        pindel -T !{task.cpus} -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}!{my_bed}
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
    
    if (params.bed) {
        INCLUDE_REGIONS(BAM_INDEX.out, params.bed)  
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
    } else {
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    }

    BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, params.bd_map_qual, params.bd_filt_qual)
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, params.pd_insert, chromosomes)
    TARDIS_PREP (BAM_INDEX.out)
    TARDIS(TARDIS_PREP.out, REFERENCE_INDEX.out, params.sonic_file, params.bed)
    MERGE_PINDEL_SINGLE(PINDEL_SINGLE.out.collect())
}
