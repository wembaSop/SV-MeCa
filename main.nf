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

process METRICS {
    
    publishDir "$params.results/metrics", mode: "copy"  

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)

    output:
        path "*.metrics"

    shell:
        outfile = "${bam.baseName}.metrics"
        '''
        MEMORY="!{task.memory}"
        MEM=${MEMORY/ GB/G}
        echo "RUNNING gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O !{outfile}"
        gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O !{outfile}
        '''

}

process BREAKDANCER {
    
    publishDir "$params.results/standalone", mode: "copy"

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        tuple path(bed), path(bed_i) //regions to include

    output:
        file "*breakdancer.vcf"
        
    shell:

        prefix = "${bam.simpleName}.breakdancer"

        '''
        run_breakdancer.sh !{bam} !{ref} !{bed} !{prefix}
        '''

}

process DELLY {
    
    publishDir "$params.results/standalone", mode: "copy" 

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
        run_delly.sh !{bam} !{ref} "!{my_bed}" "!{outfile}"
        '''

}

process LUMPY {
    
    publishDir "$params.results/standalone", mode: "copy" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed

    output:
        file "*lumpy.vcf"

    shell:

        my_bed = bed ? " --exclude $bed" : "" 
        prefix = bam.simpleName
        outfile = "${bam.simpleName}.lumpy.vcf"

        '''
        run_lumpy.sh !{bam} !{ref} !{outfile} !{prefix} !{task.cpus} "!{my_bed}"          
        
        '''

}

process MANTA {
    
    publishDir "$params.results/standalone", mode: "copy"

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
        run_manta.sh !{bam} !{ref} !{outfile} "!{task.memory}" !{task.cpus} "!{my_bed}"
        '''

}
//here update parameters
process PINDEL_SINGLE {
    
    publishDir "$params.results/pindel", mode: "copy"  

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed
        each chr

    output:
        file "*.pindel.vcf"

    shell:

        my_bed = bed ? " --exclude $bed" : ""
        fileName = bam.name
        fileSimpleName = bam.simpleName
        outfile = "${fileSimpleName}.${chr}"

        '''
        # read the sample name from bam
        SAMPLE=$(samtools samples !{bam}|head -n 1|cut -f1)

        bam2cfg.pl -g !{bam} > "!{fileSimpleName}.cfg"
        INSERT="$(grep 'readgroup' "!{fileSimpleName}.cfg"|head -n 1 | cut -f 9 |cut -d ':' -f 2)"
        INSER=${INSERT%.*}
        echo -e !{fileName}"\t$INSER\t"PD_${SAMPLE} > pindel_config.txt
        pindel -N -M 3 -r false -T !{task.cpus} -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}!{my_bed}
        pindel2vcf -P !{fileSimpleName} -is 50 -e 3 -r !{ref} -R GRCH38 -d `date +'%m/%d/%Y'` -v "!{outfile}.raw.vcf"
        bgzip "!{outfile}.raw.vcf"
        tabix -p vcf "!{outfile}.raw.vcf.gz"
        bcftools norm -f !{ref} -Ov -o "!{outfile}.norm.vcf" "!{outfile}.raw.vcf.gz"
        cat "!{outfile}.norm.vcf" | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > "!{outfile}.pindel.vcf"
        '''

}

process MERGE_PINDEL_SINGLE {
    
    publishDir "$params.results/standalone", mode: "copy" 

    input:
        path vcfs

    output:
        file "*.pindel.vcf"

    shell:

        '''
        # compress, index, concat and sort. Then add IDs 
        for f in $(ls *.vcf); do bgzip $f;tabix -p vcf "${f}.gz";done
        bcftools concat -Ov -o "${f%%.*}.pindel.unsorted.vcf" *.vcf.gz
        bcftools sort -Ov -o "${f%%.*}.pindel.raw.vcf" "${f%%.*}.pindel.unsorted.vcf"
        cat "${f%%.*}.pindel.unsorted.vcf" | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > "${f%%.*}.pindel.raw.vcf"
        grep "^#" "${f%%.*}.pindel.raw.vcf" > "${f%%.*}.pindel.vcf"
        grep -v "^#" "${f%%.*}.pindel.raw.vcf"| awk 'BEGIN{OFS="\t"};{$3="PD_"NR; print $0}'>> "${f%%.*}.pindel.vcf"

        '''

}
//To Delete
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
    
    publishDir "$params.results/standalone", mode: "copy" 

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
        run_tardis.sh !{bam} !{ref} !{sonic} !{outfile} "!{my_bed}"

        '''

}

process SURVIVOR_MERGE {
    
    publishDir "$params.results/merge", mode: "copy" 

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
        # Edit DUP to INS and ignore other sv than DEL DUP & INS
        for f in $(ls *.vcf);do edit_svtype.py $f "${f%.*}.edit.vcf";awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' "${f%.*}.edit.vcf" > "${f%.*}.edit.sort.vcf";done

        for f in $(ls *.edit.sort.vcf);do echo $f >> files.txt;done
        SURVIVOR merge files.txt 0.9 1 1 0 0 50 merge.new.vcf &> survivor.log
        cat merge.new.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.survivor.raw.vcf"
        bgzip "!{outfile}.survivor.raw.vcf"
        tabix -p vcf "!{outfile}.survivor.raw.vcf.gz"
        bcftools view -i'SVLEN<=-50 | SVLEN>=50' -Ov -o "!{outfile}.survivor.vcf" "!{outfile}.survivor.raw.vcf.gz"

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
    METRICS (BAM_INDEX.out, REFERENCE_INDEX.out)
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    if (params.pd_multi){
        chromosomes = Channel.of(1..22,"X")
        //chromosomes = Channel.of(1..22).map{"chr${it}"}
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
        MERGE_PINDEL_SINGLE(PINDEL_SINGLE.out.collect())
    } else{
        chromosomes = Channel.of("ALL")
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
    }
    
    //TARDIS_PREP (BAM_INDEX.out)
    TARDIS(BAM_INDEX.out, REFERENCE_INDEX.out, params.sonic_file, params.bed)
    SURVIVOR_MERGE(BREAKDANCER.out, DELLY.out, LUMPY.out, MANTA.out, MERGE_PINDEL_SINGLE.out, TARDIS.out)
}
