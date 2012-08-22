# Caution with orders, orders!!!!!
# Tian R., Aug. 17, 2012
cut -f1 CD133*H2AZ*dst > list
echo "gene_ID" > gene
cat gene list > headered.list

for histone in H2AZ H3K27me1 H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me1 H3K9me3 H4K20me1

do
	echo "gene	"$histone"_pro	"$histone"_trs" > header
	cat header *$histone*dst | cut -f2,3 > fea
	paste headered.list fea > temp
	cat temp > headered.list
done
rm gene
rm list
rm header
rm fea
rm temp
mv headered.list HSC.9fea.tab
	
