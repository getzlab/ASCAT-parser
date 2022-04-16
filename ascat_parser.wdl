task ascat_parser_task_1 {
    
    File caveman_tumor
    File copynumber_tumor
    File copynumber_normal
    String ID
    Int DEPTH=100
    Float HET_DENSITY_MIN=2.5e-7 
    Float ram_gb=15
    Int local_disk_gb=200
    Int num_preemptions=1

    #**Define additional inputs here**

    command {
        set -euo pipefail
        echo "$PWD" 
        ls /src
        #**Command goes here**
	python /src/ascatparser.py ${caveman_tumor} ${copynumber_tumor} ${copynumber_normal} ${DEPTH} ${HET_DENSITY_MIN} ${ID} .
    }

    output {
        #** Define outputs here**
	File ascat_allelic_capseg_file = "${ID}.ascat.acs.tsv"
        File ascat_seg_file = "${ID}.ascat.seg"
        File ascat_asc_plot = "${ID}.ascat.acs.png"
    }

    runtime {
        docker : "docker.io/chipstewart/ascat-parser_task_1:3"
        memory: "${if defined(ram_gb) then ram_gb else '4'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '20'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '1'}"
    }

    meta {
        author : "Xavi Loinaz, Alex Barbera, Chip Stewart"
        email : "xloinaz@broadinstitute.org, stewart@broadinstitute.org"
    }
}

workflow ascat_parser {

    call ascat_parser_task_1

}

