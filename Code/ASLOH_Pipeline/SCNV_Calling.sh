#!/usr/bin/bash



CNVFACETS=#cnv_facets.R path
NormalBam=#Normal Bam File
TumorBam=#Tumor Bam File
TMP=#Temporal Directory
Reference=#ucsc hg19
Output=#OutputName
Targetbed=#Target bed
dbsnp=#dbsnp_138.hg19.vcf.gz




$CNVFACETS \
-t $TumorBam \
-n $NormalBam \
-vcf $dbsnp \
--targets $Targetbed \
-g hg19 \
-o SV_$Output


