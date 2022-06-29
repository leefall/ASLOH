#!/usr/bin/bash

GATK4=#GATK4 Path
NormalBam=#Normal Bam File
TMP=#Temporal Directory
Reference=#ucsc hg19
Output=#OutputName
Targetbed=#Target bed


$GATK4 --java-options "-Xmx8G" HaplotypeCaller \
--tmp-dir=$TMP \
-R $Reference \
-I $NormalBam \
-L $Targetbed \
-O $Output.vcf.gz 


$GATK4 --java-options "-Xmx8G" SelectVariants \
--tmp-dir=$TMP \
-V $Output.vcf.gz \
-select-type SNP \
-O SNV_$Output.vcf.gz 


$GATK4 --java-options "-Xmx8G" SelectVariants \
--tmp-dir=$TMP \
-V $Output.vcf.gz \
-select-type INDEL \
-O INDEL_$Output.vcf.gz 


$GATK4 --java-options "-Xmx8G" VariantFiltration \
--tmp-dir=$TMP \
-V SNV_$Output.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "DP < 10.0" --filter-name "LowDP10" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--cluster-size 3 \
--cluster-window-size 50 \
-O Filtered_SNV_$Output.vcf.gz 


$GATK4 --java-options "-Xmx8G" VariantFiltration \
--tmp-dir=$TMP\
-V INDEL_$Output.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "DP < 10.0" --filter-name "LowDP10" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--cluster-size 3 \
--cluster-window-size 50 \
-O Filtered_INDEL_$Output.vcf.gz 













