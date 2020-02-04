# make Manhattan plots and QQplots

# pwd=/Users/leichen/Desktop/Lab/Finmetseq_paper/2-Association_test/data

SCR_plotting=/Users/leichen/Desktop/Lab/Finmetseq_paper/2-Association_test/scripts/general/6a-make_Mplot_QQplot.R


#Rscript $SCR_plotting /Users/leichen/Desktop/Lab/Finmetseq_paper/2-Association_test/data/gs/epacts_p80qt/finmetseq.ln_XXL_VLDL_TG_combined_rn.mac3.assoc.epacts

cd gs/epacts_p80qt

for f in `ls *mac3.assoc.epacts`
do 
Rscript $SCR_plotting $f
done

#mkdir 0-summary_plots
mv *pdf 0-summary_plots


cd ../../lumpy/epacts_p80qt
for f in `ls */*mac3.assoc.epacts`
do 
Rscript $SCR_plotting $f
done

#mkdir 0-summary_plots
mv */*pdf 0-summary_plots


cd ../../cnvnator/epacts_p80qt
for f in `ls */*mac3.assoc.epacts`
do 
Rscript $SCR_plotting $f
done

mkdir 0-summary_plots
mv */*pdf 0-summary_plots

# inflation:
# glu_fast 
# ln_AcAce
# S_krea
# S_ldl_c (telomere variants?)





