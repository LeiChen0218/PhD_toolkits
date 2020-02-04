x="HG03006"
sdir=/gscmnt/gc2802/halllab/aregier/jira/BIO-2978/premerge_final_outputs/pre_merge/Pre_Merge_SV_v2/b4fae28d-ee3f-4614-8688-8343eb862a75/call-Pre_Merge_SV_Per_Sample/shard-1144/per_sample.Pre_Merge_SV_Per_Sample/fa9ee881-1221-4481-bf55-db369bbd4603/call-CNVnator_Histogram/cnvnator.out/HG03006.final.cram.hist.root
F_DIR=/gscmnt/gc2802/halllab/wangciyang/projects/data #where you have the window file generated
WINDOW_SIZE=1000 #change to your own window size
scr_cnvnator=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/3-Candidate_analysis/scripts/cnvnator_b38.full_name.sh

LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
bsub -a 'docker(halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2)' \
-oo ${F_DIR}/cn/logs/$x.cnvnator.log \
-q ccdg -n 5 -M 50000000 \
-R "rusage[gtmp=10, mem=50000] select[mem>50000]" \
-g /lchen/cnvnator \
/bin/bash ${scr_cnvnator} $x $sdir ${F_DIR}/region.${WINDOW_SIZE}bp.windows.list ${F_DIR}
