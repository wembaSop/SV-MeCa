//support alignment
//support CRAM input
// support VCF inputs

process REFERENCE_INDEX {
    
    publishDir "$params.results/reference_index", enabled:false 

    input:
        path ref

    output:
        tuple path(ref), path("*.fai")

    shell:

        '''
        samtools faidx !{ref}
        '''
    stub:
        """
        samtools version > test.fa.fai
        """

}

process INCLUDE_REGIONS {

    input:
        tuple path(bam), path(bai) 
        path bed

    output:
        tuple path("include.bed.gz"), path("include.bed.gz.tbi")

    shell:

        '''
        samtools view -H !{bam}| grep "^@SQ"| cut -f 2-3|sed -e 's/SN://g' -e 's/LN://g'|awk '{print $1"\t"0"\t"$2}'> full.bed
        bedtools subtract -a full.bed -b !{bed} > include.bed
        bgzip include.bed
        tabix -p bed include.bed.gz

        '''
    stub:
        """
        samtools version > full.bed
        bedtools --version > include.bed
        bgzip --version > include.bed.gz
        tabix --version > include.bed.gz.tbi

        """

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
    stub:
        """
        samtools version > test.bam.bai
        """
}

process METRICS {
    
    publishDir "$params.results/metrics", mode: "copy"  

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)

    output:
        path "*.metrics"

    shell:
        outfile = "${bam.simpleName}.metrics"
        '''
        MEMORY="!{task.memory}"
        MEM=${MEMORY/ GB/G}
        echo "RUNNING gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O !{outfile}"
        gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O !{outfile}
        '''

    stub:

        outfile = "${bam.simpleName}.metrics"
    
        """
        MEMORY="${task.memory}"
        MEM=\${MEMORY/ GB/G}
        echo "RUNNING gatk --java-options "-Xmx\$MEM" CollectWgsMetrics -I ${bam} -R ${ref} -O ${outfile}" > $outfile
        gatk CollectWgsMetrics --version >> $outfile
        """

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
    stub:

        prefix = "${bam.simpleName}.breakdancer"

        """
        touch "${prefix}.cfg"
        touch "${prefix}.ctx"
        breakdancertovcf.py -h > "${prefix}.raw.vcf"
        grep -V > "${prefix}.unsorted.vcf"
        awk -V >> "${prefix}.unsorted.vcf"
        bgzip --version > "${prefix}.vcf.gz"
        tabix --version > "${prefix}.vcf.gz.zbi"
        bcftools -v > "${prefix}.vcf"
        sed --version >> "${prefix}.vcf"
        """

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
    stub:

        my_bed = bed ? " -x $bed" : "" 
        outfile = "${bam.simpleName}.delly.vcf"

        """
        delly --help > $outfile
        bcftools -v >> ${outfile}
        sed --version >> ${outfile}
        """

}

process INSURVEYOR {
    
    publishDir "$params.results/standalone", mode: "copy" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        tuple path(bed), path(bed_i)

    output:
        file "*insurveyor.vcf"

    shell:
        outfile = "${bam.simpleName}.insurveyor.vcf"

        '''
        SAMPLE=$(samtools samples !{bam}|head -n 1|cut -f1)
        insurveyor.py --threads !{task.cpus} --samplename "IS_${SAMPLE}" --min-insertion-size 50 !{bam} $PWD !{ref}
        bcftools view -R !{bed} -Ov -o !{outfile} out.vcf.gz
        '''
    stub:
        outfile = "${bam.simpleName}.insurveyor.vcf"
        """
        insurveyor.py -h &> $outfile
        bcftools -v >> $outfile
        sed --version >> $outfile
        """

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
    stub:

        my_bed = bed ? " --exclude $bed" : "" 
        prefix = bam.simpleName
        outfile = "${bam.simpleName}.lumpy.vcf"

        """
        smoove call --help > file.lumpy.vcf

        """

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
    stub:

        """
        touch file.manta.vcf
        """

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
    stub:
        my_bed = bed ? " --exclude $bed" : ""
        fileName = bam.name
        fileSimpleName = bam.simpleName
        outfile = "${fileSimpleName}.${chr}"
        """
        touch "${outfile}.pindel.vcf"
        """

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
    stub:
        """
        touch file.pindel.vcf
        """

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
    stub:
        """
        touch file.markdup.bam
        touch file.markdup.bam.bai
        """

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

        my_bed = bed ? " --gaps $bed" : " " 
        outfile = "${bam.simpleName}"

        '''
        run_tardis.sh !{bam} !{ref} !{sonic} !{outfile} "!{my_bed}"

        '''

    stub:

        """
        touch file.tardis.vcf
        """

}

process SURVIVOR_MERGE {
    
    publishDir "$params.results/merge", mode: "copy" 

    input:
        path breakdancer
        path delly
        path insurveyor
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
        bcftools view -i 'SVLEN<=-50 | SVLEN>=50' -Ov -o "!{outfile}.survivor.vcf" "!{outfile}.survivor.raw.vcf.gz"

        '''

    stub:
    
        """
        touch file.survivor.vcf
        """

}



workflow {
    
    REFERENCE_INDEX (params.reference)
    BAM_INDEX (params.input)

    if (params.bed) {
        INCLUDE_REGIONS(BAM_INDEX.out, params.bed)  
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
        INSURVEYOR (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
    } else {
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
        INSURVEYOR (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    }
    
    METRICS (BAM_INDEX.out, REFERENCE_INDEX.out)
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed)

    if (params.pd_multi){
        //chromosomes = Channel.of(1..22,"X")
        chromosomes = Channel.of(1..22).map{"${it}"}
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
        MERGE_PINDEL_SINGLE(PINDEL_SINGLE.out.collect())
    } else{
        chromosomes = Channel.of("ALL")
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, params.bed, chromosomes)
    }

    if (params.enable_markdup){
        TARDIS_PREP (BAM_INDEX.out)
        TARDIS(TARDIS_PREP.out, REFERENCE_INDEX.out, params.sonic, params.bed)
    } else {
        TARDIS(BAM_INDEX.out, REFERENCE_INDEX.out, params.sonic, params.bed)
    }
    
    SURVIVOR_MERGE(BREAKDANCER.out, DELLY.out, INSURVEYOR.out, LUMPY.out, MANTA.out, MERGE_PINDEL_SINGLE.out, TARDIS.out)
}
