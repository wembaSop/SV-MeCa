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
        tuple path("include.bed.gz"), path("include.bed.gz.tbi")

    shell:

        '''
        samtools view -H !{bam}| grep "^@SQ"| awk -F '[\t:]' '{print $3"\t"0"\t"$5}'> full.bed
        bedtools subtract -a full.bed -b !{bed} > include.bed
        bgzip include.bed
        tabix -p bed include.bed.gz

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
    
    publishDir "$params.results/standalone" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        tuple path(bed), path(bed_i) //regions to include

    output:
        file "*breakdancer.vcf"
        
    shell:

        prefix = "${bam.simpleName}.breakdancer"

        '''
        bam2cfg.pl -g !{bam} > "!{prefix}.cfg"
        breakdancer-max -h "!{prefix}.cfg" > "!{prefix}.ctx"
        breakdancertovcf.py -o "!{prefix}.raw.vcf" !{ref} "!{prefix}.ctx"
        grep "^#" "!{prefix}.raw.vcf" > "!{prefix}.unsorted.vcf"
        grep -v "^#" "!{prefix}.raw.vcf"| awk 'BEGIN{OFS="\t"};{$3=NR; print $0}'>> "!{prefix}.unsorted.vcf"
        cat "!{prefix}.unsorted.vcf" | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > "!{prefix}.sorted.vcf"
        bgzip "!{prefix}.sorted.vcf"
        tabix -p vcf "!{prefix}.sorted.vcf.gz"
        bcftools view -R !{bed} -Ov -o "!{prefix}.vcf" "!{prefix}.sorted.vcf.gz"
        '''

}

process DELLY {
    
    publishDir "$params.results/standalone" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed 

    output:
        file "*delly.vcf"

    shell:

        my_bed = bed ? " -x $bed" : "" 
        outfile = "${bam.simpleName}.delly.vcf"

        '''
        delly call!{my_bed} -t ALL -g !{ref} -o delly_output.bcf !{bam}
        bcftools view -Ov -o !{outfile} delly_output.bcf 
        '''

}

process LUMPY {
    
    publishDir "$params.results/standalone" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed

    output:
        file "*lumpy.vcf"

    shell:

        my_bed = bed ? " --exclude $bed" : "" 
        prefix = bam.simpleName
        outfile = "${bam.simpleName}.lumpy.vcf.gz"

        '''
        smoove call -x --genotype --name !{prefix} --outdir . -f !{ref} --processes !{task.cpus}!{my_bed} !{bam}
        mv *genotyped.vcf.gz !{outfile}
        bgzip -d !{outfile} 
        
        '''

}

process MANTA {
    
    publishDir "$params.results/standalone"

    input:
        tuple path(bam), path(bam_index)
        tuple path(ref), path(ref_index)
        tuple path(bed), path(bed_i) //regions to include

    output:
        file "*manta.vcf"

    shell:

        my_bed = bed ? " --callRegions $bed --exome" : ""  
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
        each chr

    output:
        file "*.vcf"

    shell:

        my_bed = bed ? " --exclude $bed" : ""
        fileName = bam.name
        fileSimpleName = bam.simpleName
        outfile = "${fileSimpleName}.pindel.${chr}.vcf"

        '''
        bam2cfg.pl -g !{bam} > "!{fileSimpleName}.cfg"
        INSERT="$(grep 'readgroup' "!{fileSimpleName}.cfg"|head -n 1 | cut -f 9 |cut -d ':' -f 2)"
        INSER=${INSERT%.*}
        echo -e !{fileName}"\t$INSER\t"!{fileSimpleName} > pindel_config.txt
        pindel -N -M 3 -r false -T !{task.cpus} -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}!{my_bed}
        pindel2vcf -P !{fileSimpleName} -is 50 -e 3 -r !{ref} -R GRCH37 -d `date +'%m/%d/%Y'` -v !{outfile}

        '''

}

process MERGE_PINDEL_SINGLE {
    
    publishDir "$params.results/standalone" 

    input:
        path vcfs

    output:
        file "*.pindel.vcf"

    shell:

        '''
        for f in $(ls *.vcf); do bgzip $f;tabix -p vcf "${f}.gz";done
        bcftools concat -Ov -o "${f%%.*}.pindel.unsorted.vcf" *.vcf.gz
        bcftools sort -Ov -o "${f%%.*}.pindel.raw.vcf" "${f%%.*}.pindel.unsorted.vcf"
        grep "^#" "${f%%.*}.pindel.raw.vcf" > "${f%%.*}.pindel.vcf"
        grep -v "^#" "${f%%.*}.pindel.raw.vcf"| awk 'BEGIN{OFS="\t"};{$3=NR; print $0}'>> "${f%%.*}.pindel.vcf"
        '''

}

process TARDIS_PREP {
    
    //publishDir "$params.results" 

    input:
        tuple path(bam),path(bam_index)

    output:
        tuple path("*.markdup.bam"), path("*.markdup.bam.bai")

    shell:

        outfile = "${bam.simpleName}.markdup.bam"

        '''
        mkdir tmp_sam
        sambamba markdup -r -t !{task.cpus} --tmpdir=tmp_sam !{bam} !{outfile}
        samtools index -@!{task.cpus} !{outfile}

        '''

}

process TARDIS {
    
    publishDir "$params.results/standalone" 

    input:
        tuple path(bam), path(bam_index)
        tuple path(ref),path(ref_index)
        path sonic
        path bed 

    output:
        file "*.tardis.vcf"

    shell:

        my_bed = bed ? " --gaps $bed" : "" 
        outfile = "${bam.baseName}"

        '''
        tardis -i !{bam} --ref !{ref} --sonic !{sonic} --out !{outfile}!{my_bed}
        cat "!{outfile}.vcf" | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.tardis.vcf"

        '''

}

process SURVIVOR_MERGE {
    
    publishDir "$params.results/merge" 

    input:
        path breakdancer
        path delly
        path lumpy
        path manta
        path pindel
        path tardis

    output:
        file "*.survivor.vcf"

    shell:

        
        outfile = "${delly.simpleName}"

        '''
        for f in $(ls *.vcf);do edit_svtype.py $f "${f%.*}.edit.vcf";done
        for f in $(ls *.edit.vcf);do echo $f >> files.txt;done
        SURVIVOR merge files.txt 0.9 1 1 0 0 50 merge.new.vcf &> survivor.log
        cat merge.new.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.survivor.raw.vcf"
        bgzip "!{outfile}.survivor.raw.vcf"
        tabix -p vcf "!{outfile}.survivor.raw.vcf.gz"
        bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -Ov -o "!{outfile}.survivor.vcf" "!{outfile}.survivor.raw.vcf.gz"

        '''

}



workflow {
    
    REFERENCE_INDEX (params.reference)
    BAM_INDEX (params.input)
    if (params.bed) {
        INCLUDE_REGIONS(BAM_INDEX.out, params.bed)  
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
    } else {
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    }
    
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    if (params.pd_multi){
        chromosomes = Channel.of(1..22)
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
        MERGE_PINDEL_SINGLE(PINDEL_SINGLE.out.collect())
    } else{
        chromosomes = Channel.of("ALL")
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
    }
    
    TARDIS_PREP (BAM_INDEX.out)
    TARDIS(TARDIS_PREP.out, REFERENCE_INDEX.out, params.sonic_file, params.bed)
    SURVIVOR_MERGE(BREAKDANCER.out, DELLY.out, LUMPY.out, MANTA.out, MERGE_PINDEL_SINGLE.out, TARDIS.out)
}
