task ascat_parser_task_1 {
    
    Float ram_gb = 15
    Int local_disk_gb = 200
    Int num_preemptions = 1

    #**input files from ascatNGS**
    File ascat_caveman_file
    File ascat_copynumber_file
    File ascat_copynumber_normal_file
    String ID


    command {
        set -euo pipefail
		echo $pwd
        echo "ascat_caveman_file: ${ascat_caveman_file}"
        echo "ascat_copynumber_file: ${ascat_copynumber_file}"
        echo "ascat_copynumber_normal_file: ${ascat_copynumber_normal_file}"
        echo "ID: ${ID}"

    	python ./ascatparser.py ${ascat_caveman_file} ${ascat_copynumber_file} ${ascat_copynumber_normal_file} ${ID}.ascat.allelic_capseg.tsv ${ID}.ascat.seg
 
    }

    output {
        File ascat_allelic_capseg_tsv = "${ID}.ascat.allelic_capseg.tsv" 
        File ascat_seg = "${ID}.ascat.seg" 
    }

    runtime {
        docker : "gcr.io/xloinaz/ascat_parser_task_1:1"
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
    
    
