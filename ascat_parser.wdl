task ascat_parser_task_1 {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    #**Define additional inputs here**

    command {
        set -euo pipefail

        #**Command goes here**
    }

    output {
        #** Define outputs here**
    }

    runtime {
        docker : "<namespace>/ascat_parser_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : ""
        email : "xloinaz@broadinstitute.org"
    }
}

workflow ascat_parser {

    call ascat_parser_task_1 {
        input: #**Define call inputs for ascat_parser_task_1 here**
    }

    output {
        #**Define workflow outputs here. If defined, these will be the only
        #  outputs available in the Method Configuration**
    }
}
