INPUT_VCF=$1
SNP_VCF=$2
evalDir=$3

JAVA_HOME=/gapp/x64linux/opt/java/jdk/jdk1.8.0_60
JAVA=${JAVA_HOME}/bin/java

EXPERIMENTAL_GATK=/gscmnt/gc2802/halllab/dlarson/jira/BIO-2633/vendor/local/jars/
GATK=${EXPERIMENTAL_GATK}/GenomeAnalysisTK-3.8-1-idas-experimental-4131062-2018.01.27.jar
QUEUE=${EXPERIMENTAL_GATK}/Queue-3.8-1-idas-experimental-07c96d8-2018.02.17.jar

GENOMESTRIP=/gscmnt/gc2802/halllab/ccdg_resources/pkg/GenomeSTRiP_2.00.1774
SV_DIR=${GENOMESTRIP}/svtoolkit
SVTOOLKIT=${SV_DIR}/lib/SVToolkit.jar
CLASSPATH="${SVTOOLKIT}:${GATK}:${QUEUE}"

METADATA_DIR=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/GenomeSTRiP/Homo_sapiens_assembly38
REF=${METADATA_DIR}/Homo_sapiens_assembly38.fasta

CONFIG=${SV_DIR}/conf/genstrip_parameters.txt
PLOIDY_MAP=${METADATA_DIR}/Homo_sapiens_assembly38.ploidymap.txt

$JAVA -cp $CLASSPATH -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A TagVariants \
    -R ${REF} \
    -vcf ${INPUT_VCF} \
    -configFile ${CONFIG} \
    -tagFile ${SNP_VCF} \
    -ploidyMapFile ${PLOIDY_MAP} \
    -writeReport true \
    -writeSummary true \
    -reportDirectory ${evalDir} \
    || exit 1