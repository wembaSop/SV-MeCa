
process GENOME_INDEX {
    
    publishDir "$params.results/genome_index" 
    input:
        path genome
    output:
        path "${genome}.fai"
    shell:
        '''
        samtools faidx !{genome}
        '''

}

process BAM_INDEX {
    
    publishDir "$params.results/bam_index" 
    input:
        path bam
    output:
        path "${bam}.bai"
    shell:
        '''
        samtools index !{bam}
        '''

}

process BREAKDANCER {
    
    publishDir "$params.results/breakdancer" 
    input:
        path bam
        path bam_index
        path ref
    output:
        file "*.ctx"
        
    shell:
        prefix = bam.baseName
        '''
        bam2cfg.pl -q 20 -g !{bam} > "!{prefix}.cfg"
        breakdancer-max -q 20 -y 20 -h "!{prefix}.cfg" > "!{prefix}.ctx"
        python breakdancertovcf.py -o "!{prefix}.vcf" !{ref} "!{prefix}.ctx"
        
        '''

}


process MANTA {
    
    publishDir "$params.results/manta" 
    input:
        path bam
        path bam_index
        path ref 
        path ref_index
    output:
        file "*.vcf.gz"
    shell:
        '''
        configManta.py --bam !{bam} --referenceFasta !{ref} --runDir ./
        ./runWorkflow.py
        cp results/variants/diploidSV.vcf.gz .
        '''

}
process LUMPY {
    
    publishDir "$params.results/lumpy" 
    input:
        path bam
        path bam_index
        path bed
        path ref 
        path ref_index
    output:
        file "*.vcf.gz"
    shell:
        '''

        smoove call -x --genotype --name hg002_smoove --outdir . -f !{ref} --processes 4 --exclude !{bed} !{bam}
        
        '''

}
process DELLY {
    
    publishDir "$params.results/delly" 
    input:
        path bam
        path bam_index
        path ref 
        path ref_index
    output:
        file "*.vcf.gz"
    shell:
        '''
        delly call -t ALL -g !{ref} -o delly_output.bcf !{bam}
        bcftools view -Oz -o delly_output.vcf.gz delly_output.bcf 
        '''

}

process PINDEL_SINGLE {
    
    publishDir "$params.results/pindel" 
    input:
        path bam
        path bam_index
        path ref 
        path ref_index
        each chr
    output:
        file "*.vcf"
    shell:
    fileName = bam.name
    fileSimpleName = bam.simpleName
    outfile = "${fileSimpleName}.${chr}.vcf"
        '''
        echo -e !{fileName}'\t'350'\t'!{fileSimpleName} > pindel_config.txt
        pindel -T 12 -x 5 -f !{ref} -i pindel_config.txt -c !{chr} -o !{fileSimpleName}
        pindel2vcf -P !{fileSimpleName} -r !{ref} -R GRCH37 -d `date +'%m/%d/%Y'` -v !{outfile}

        '''

}

process SAMBAMBA {
    
    publishDir "$params.results/sambamba" 
    input:
        path bam
        path bam_index

    output:
        tuple path("*.markdup.bam"), path("*.markdup.bam.bai")
    shell:
        outfile = "${bam.baseName}.markdup.bam"
        '''
        mkdir tmp_sam
        sambamba markdup -r --tmpdir=tmp_sam !{bam} !{outfile}
        sambamba index !{outfile}

        '''

}

process TARDIS {
    
    publishDir "$params.results/tardis" 
    input:
        tuple path(bam), path(bam_index)
        path ref 
        path ref_index
        path sonic
    output:
        file "*.vcf"
    shell:
        outfile = "${bam.baseName}"
        '''
        tardis -i !{bam} --ref !{ref} --sonic !{sonic} --out !{outfile}
        cat "!{outfile}.vcf" | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > "!{outfile}.sorted.vcf"

        '''

}
workflow {
    CHR = [1,2,321,22]
    //GENOME_INDEX (params.reference)
    BAM_INDEX (params.input)
    //MANTA (params.input, BAM_INDEX.out, params.reference, GENOME_INDEX.out)
    //LUMPY (params.input, BAM_INDEX.out,params.exclude_bed, params.reference, GENOME_INDEX.out)
    //DELLY (params.input, BAM_INDEX.out, params.reference, GENOME_INDEX.out)
    //PINDEL_SINGLE (params.input, BAM_INDEX.out, params.reference, GENOME_INDEX.out, CHR)
    //SAMBAMBA (params.input, BAM_INDEX.out)
    //TARDIS(SAMBAMBA.out, params.reference, GENOME_INDEX.out, params.sonic_file)
    BREAKDANCER(params.input, BAM_INDEX.out)
}
