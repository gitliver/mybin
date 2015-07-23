#!/bin/bash

input_bam=${1} 
output_bam=${2}
metrics_file=${3}
java_memory=${4}
picard=${5}
# tmp=/path/to/tmp/dir

### VALIDATION_STRINGENCY=LENIENT, that will tell picard to show any error it sees but to continue with the processing
# java -Xmx${java_memory}G -Djava.io.tmpdir=$tmp -jar $picard/MarkDuplicates.jar \
java -Xmx${java_memory}G -jar $picard/MarkDuplicates.jar \
    INPUT=$input_bam \
    OUTPUT=$output_bam \
    METRICS_FILE=$metrics_file \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT
