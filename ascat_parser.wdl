task ascat_parser_task_1 {
    
    File caveman_tumor
    File copynumber_tumor
    File copynumber_normal
    String ID
 
    Float ram_gb=15
    Int local_disk_gb=200
    Int num_preemptions=1

    #**Define additional inputs here**

    command {
        set -euo pipefail
        echo "$PWD" 
        ls /src
        #**Command goes here**
		python /src/ascatparser.py ${caveman_tumor} ${copynumber_tumor} ${copynumber_normal} ${ID}.ascat_allelic_capseg.tsv ${ID}.ascat.seg 0 
    }

    output {
        #** Define outputs here**
		File ascat_allelic_capseg_file = "${ID}.ascat_allelic_capseg.tsv"
        File ascat_gistic_input_file = "${ID}.ascat.seg"
    }

    runtime {
        docker : "gcr.io/broad-getzlab-wolf-wgs-hg38/xloinaz/ascat-parser_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Xavi Loinaz"
        email : "xloinaz@broadinstitute.org"
    }
}

workflow ascat_parser {

    call ascat_parser_task_1

}

