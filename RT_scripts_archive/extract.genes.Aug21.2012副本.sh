# Tian R. Aug 21, 2012
# extract genes according to their expression changes in two cell types
# HSC ERY

for file in `ls *exprs.tab`

do

	head -n1 HSC_9fea_30k_NM_Aug17_2012.tab  > $file".9fea" 
						#stay_silen_9fea.TAB

	#silen2silen.exprs.tab 
	cat $file | sed 's/"//g' | cut -f1 | while IFS=\t read gene 
		do
			echo $gene
			grep -w $gene HSC_9fea_30k_NM_Aug17_2012.tab >> $file".9fea"
		done
done


