# the master script which runs other scripts to generate the plot

ROOT=/gscmnt/gc2802/halllab/lchen/cn_plot_demo  # your own project directory 
cd $ROOT
# mkdir data scripts logs
# transfer 1-cn_window_1kg.sh, 2-check_output_1kg.sh, and 3-paste_cn_files.sh to ${ROOT}/scripts/
DIR=${ROOT}/data  # your working directory

cd $DIR

# here's an example for plotting the deletion chr11_26969501_27010500
var="chr11_26969501_27010500" # specfy your own variant of interests
chr=11
pos=26966201
end=27220000

window_size=1000
window_number=200

# 1.Generate windows for the variant and its flanking regions (window size: 1kb, 200 windows on each side)
window_scr=$ROOT/scripts/1-generate_windows.sh
bash $window_scr $chr $pos	$end $window_size $window_number $DIR 2> ${ROOT}/logs/1-generate_windows.log

# 2.Extract copy number data for each sample (20 examples)
# 	cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/3-Candidate_analysis/data/var_plots/chr11_26969501_27010500/data 
#	awk '$2==1' genotype.discrete.t.list | head > /gscmnt/gc2802/halllab/lchen/cn_plot_demo/data/exp.20samples.list
#	awk '$2==4' genotype.discrete.t.list | head >> /gscmnt/gc2802/halllab/lchen/cn_plot_demo/data/exp.20samples.list
#	cd $DIR
# 	sv_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.sample.all.table
# 	cut -f 1,5 ${sv_dir} | zjoin -a stdin -b exp.20samples.list -wa > exp.20samples.dir.list
cn_scr=$ROOT/scripts/2-extract_cn.sh
bash $cn_scr exp.20samples.dir.list $DIR/flanking_${window_size} $window_size

# 3.check the results to see if every sample finished
check_scr=$ROOT/scripts/3-check_output.sh
bash $check_scr $DIR/flanking_${window_size} exp.20samples.dir.list 2> ${ROOT}/logs/3-check_output.log

# 4.paste the results from all samples and generate the matrix
paste_scr=$ROOT/scripts/4-paste_cn_files.sh
bash $paste_scr $DIR/flanking_${window_size} exp.20samples.dir.list 2> ${ROOT}/logs/4-paste_cn_files.log
awk '{print NF}' $DIR/flanking_${window_size}/all_sample.cn  | sort | uniq -c # check the dimention of the matrix

########
# 5. pull the all_sample.cn file down and make cn plots on local laptop
# Some R packages might be required


