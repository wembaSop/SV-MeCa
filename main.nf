//support alignment
//support CRAM input
// support VCF inputs

PUBLISH = "$params.results"

process REFERENCE_INDEX {
    
    publishDir "$PUBLISH/reference_index", enabled:false 

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
    
    //publishDir "$PUBLISH/bam_index" 

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
    
    publishDir "$PUBLISH/metrics", mode: "copy"  

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path cfg

    output:
        path "*.metrics", emit: metrics
        path "*.stats", emit: stats

    shell:
        outfile = "${bam.simpleName}"
        '''
        MEMORY="!{task.memory}"
        MEM=${MEMORY/ GB/G}
        echo "RUNNING gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O !{outfile}.metrics"
        gatk --java-options "-Xmx$MEM" CollectWgsMetrics -I !{bam} -R !{ref} -O "!{outfile}.metrics"
        echo "readlen=$(grep 'readgroup' !{cfg} | head -n 1 | cut -f 4 |cut -d ':' -f 2)" >> "!{outfile}.metrics.stats"
        echo "coverage=$(head -n 8 !{outfile}.metrics |tail -n1 |cut -f 2)" >> "!{outfile}.metrics.stats"
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
    
    publishDir "$PUBLISH/standalone", mode: "copy"

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        tuple path(bed), path(bed_i) //regions to include

    output:
        path "*breakdancer.vcf", emit: vcf
        path "*breakdancer.cfg", emit: cfg
        
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
    
    publishDir "$PUBLISH/standalone", mode: "copy" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed 

    output:
        path "*delly.vcf"

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
    
    publishDir "$PUBLISH/standalone", mode: "copy" 

    input:
        tuple path(bam), path(bam_index)
        tuple path(ref), path(ref_index)
        tuple path(bed), path(bed_i)

    output:
        path "*insurveyor.vcf"

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
    
    publishDir "$PUBLISH/standalone", mode: "copy" 

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed

    output:
        path "*lumpy.vcf"

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
    
    publishDir "$PUBLISH/standalone", mode: "copy"

    input:
        tuple path(bam), path(bam_index)
        tuple path(ref), path(ref_index)
        tuple path(bed), path(bed_i) //regions to include

    output:
        path "*manta.vcf"

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
    
    publishDir "$PUBLISH/pindel", mode: "copy"  

    input:
        tuple path(bam),path(bam_index)
        tuple path(ref),path(ref_index)
        path bed
        each chr

    output:
        path "*.pindel.vcf"

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
        pindel -N -M 3 -r false -w 1 -T !{task.cpus} -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}!{my_bed}
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
    
    publishDir "$PUBLISH/standalone", mode: "copy" 

    input:
        path vcfs

    output:
        path "*.pindel.vcf"

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
    
    //publishDir "$PUBLISH" 

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
    
    publishDir "$PUBLISH/standalone", mode: "copy" 

    input:
        tuple path(bam), path(bam_index)
        tuple path(ref),path(ref_index)
        path sonic
        path bed 

    output:
        path "*.tardis.vcf"

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
    
    publishDir "$PUBLISH/merge", mode: "copy" 

    input:
        path breakdancer
        path delly
        path insurveyor
        path lumpy
        path manta
        path pindel
        path tardis
        path stats

    output:
        path "*.survivor.vcf", emit: survivor
        path "*.breakdancer.edit.sort.vcf", emit: breakdancer
        path "*.delly.edit.sort.vcf", emit: delly
        path "*.insurveyor.edit.sort.vcf", emit: insurveyor
        path "*.lumpy.edit.sort.vcf", emit: lumpy
        path "*.manta.edit.sort.vcf", emit: manta
        path "*.pindel.edit.sort.vcf", emit: pindel
        path "*.tardis.edit.sort.vcf", emit: tardis
        path "*.stats", emit: stats


    shell:

        
        outfile = "${delly.simpleName}"

        '''
        cp !{stats} "!{outfile}.merge.stats"
        # Edit DUP to INS and ignore other sv than DEL DUP & INS
        for f in $(ls *.vcf);do tool=${f#*.}; echo "${tool%%.*}=$(grep -cv "^#" $f)" >> "!{outfile}.merge.stats"; edit_svtype.py $f "${f%.*}.edit.vcf";awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' "${f%.*}.edit.vcf" > "${f%.*}.edit.sort.vcf";done

        for f in $(ls *.edit.sort.vcf);do echo $f >> files.txt;done
        SURVIVOR merge files.txt 0.9 1 1 0 0 50 merge.new.vcf &> survivor.log
        cat merge.new.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.survivor.raw.vcf"
        bgzip "!{outfile}.survivor.raw.vcf"
        tabix -p vcf "!{outfile}.survivor.raw.vcf.gz"
        bcftools view -i 'SVLEN<=-50 | SVLEN>=50' -Ov -o "!{outfile}.survivor.vcf" "!{outfile}.survivor.raw.vcf.gz"
        echo "survivor=$(grep -cv "^#" !{outfile}.survivor.vcf)" >> "!{outfile}.merge.stats"
        '''

    stub:
    
        """
        touch file.survivor.vcf
        """

}

process SCORING {
    
    publishDir "$PUBLISH/sv-meca", mode: "copy" 

    input:
        path breakdancer
        path delly
        path insurveyor
        path lumpy
        path manta
        path pindel
        path tardis
        path survivor
        path model_del
        path model_ins
        path stats

    output:
        path "*.svmeca.vcf"

    shell:

        outfile = "${delly.simpleName}"

        '''
        sv_meca.py -o "!{outfile}.svmeca.raw.vcf" -s !{stats} -bd !{breakdancer} -dl !{delly} -is !{insurveyor} -lp !{lumpy} -mt !{manta} -pd !{pindel} -td !{tardis} -su !{survivor} !{model_ins} !{model_del}
        awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' "!{outfile}.svmeca.raw.vcf" > "!{outfile}.svmeca.vcf"
        '''

    stub:

        """
        touch file.svmeca.vcf
        """

}


workflow sv_calling {
    take:
    bam_file
    ref_file
    bed_file
    chr_wise
    has_chr
    sonic
    td_markdup

    
    main:
    REFERENCE_INDEX (ref_file)
    BAM_INDEX (bam_file)

    if (bed_file) {
        INCLUDE_REGIONS(BAM_INDEX.out, bed_file)  
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
        INSURVEYOR (BAM_INDEX.out, REFERENCE_INDEX.out, INCLUDE_REGIONS.out)
    } else {
        MANTA (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file)
        BREAKDANCER (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file)
        INSURVEYOR (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file)
    }
    
    METRICS (BAM_INDEX.out, REFERENCE_INDEX.out, BREAKDANCER.out.cfg)
    DELLY (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file)
    LUMPY (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file)

    if (chr_wise){
        if (has_chr){
            //chromosomes = Channel.of(1..22,"X")
            chromosomes = Channel.of(1..22).map{"chr${it}"}
            PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file, chromosomes)
            MERGE_PINDEL_SINGLE (PINDEL_SINGLE.out.collect())
        } else {
            //chromosomes = Channel.of(1..22,"X")
            chromosomes = Channel.of(1..22).map{"${it}"}
            PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file, chromosomes)
            MERGE_PINDEL_SINGLE (PINDEL_SINGLE.out.collect())
        }

    } else{
        chromosomes = Channel.of("ALL")
        PINDEL_SINGLE (BAM_INDEX.out, REFERENCE_INDEX.out, bed_file, chromosomes)
    }

    if (td_markdup){
        TARDIS_PREP (BAM_INDEX.out)
        TARDIS (TARDIS_PREP.out, REFERENCE_INDEX.out, sonic, bed_file)
    } else {
        TARDIS (BAM_INDEX.out, REFERENCE_INDEX.out, sonic, bed_file)
    }

    emit:
    breakdancer = BREAKDANCER.out.vcf
    delly = DELLY.out
    insurveyor = INSURVEYOR.out
    lumpy = LUMPY.out
    manta = MANTA.out
    pindel = MERGE_PINDEL_SINGLE.out
    tardis = TARDIS.out
    stats = METRICS.out.stats

}

workflow merge_score{
    take: 
    breakdancer
    delly
    insurveyor
    lumpy
    manta
    pindel
    tardis
    model_del
    model_ins
    stats

    main:
    SURVIVOR_MERGE (breakdancer, delly, insurveyor, lumpy, manta, pindel, tardis, stats)
    SCORING (SURVIVOR_MERGE.out.breakdancer, SURVIVOR_MERGE.out.delly, SURVIVOR_MERGE.out.insurveyor, SURVIVOR_MERGE.out.lumpy, SURVIVOR_MERGE.out.manta, SURVIVOR_MERGE.out.pindel, SURVIVOR_MERGE.out.tardis, SURVIVOR_MERGE.out.survivor, model_del, model_ins, SURVIVOR_MERGE.out.stats)

    emit:
    sv_meca_vcf = SCORING.out

}

workflow{
    // Start from FASTQ files.  
    if (params.format == "fastq"){
        
        //TO DO
        println "FASTQ COMING SOON"

    // Start from BAM files. 
    } else if (params.format == "bam"){
        
        sv_calling(params.input, params.reference, params.bed, params.pd_multi, params.has_chr, params.sonic, params.enable_markdup)
        merge_score(sv_calling.out.breakdancer, sv_calling.out.delly, sv_calling.out.insurveyor, sv_calling.out.lumpy, sv_calling.out.manta, sv_calling.out.pindel, sv_calling.out.tardis, params.model_del, params.model_ins, sv_calling.out.stats)
    
    //start from VCF
    } else if (params.format == "vcf"){

       merge_score(params.breakdancer, params.delly, params.insurveyor, params.lumpy, params.manta, params.pindel, params.tardis, params.model_del, params.model_ins, params.stats)
    
    } else {

        println "Provide the parameter format via --format or in the config file"
    
    }

}
