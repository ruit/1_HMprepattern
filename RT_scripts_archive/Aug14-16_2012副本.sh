#axel -n20 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/GSE12646/GSE12646_RAW.tar

#Aug 14, 2012
#CD133-H2AZ
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005144/SRR016977/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005144/SRR016978/*.sra

#CD133-H3K27me1
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005145/SRR016979/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005145/SRR016980/*.sra

#CD133-H3K27me3
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005146/SRR016981/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005146/SRR016982/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005146/SRR016983/*.sra

#CD133-H3K36me3
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005147/SRR016984/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005147/SRR016985/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005147/SRR016986/*.sra

#CD133-H3K4me1
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005148/SRR016987/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005148/SRR016988/*.sra 

#CD133-H3K4me3
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005149/SRR016989/*.sra

#CD133-H3K9me1
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005150/SRR016990/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005150/SRR016991/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005150/SRR016992/*.sra

#CD133-H3K9me3
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005151/SRR016993/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005151/SRR016994/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005151/SRR016995/*.sra

#CD133-H4K20me1
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005152/SRR016996/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005152/SRR016997/*.sra
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005152/SRR016998/*.sra

#CD133-PolII
#axel -n20 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX005/SRX005153/SRR016999/*.sra

#
#fastq-dump

for file in `ls *sra`
do
	echo $file
	sratoolkit.2.1.16-centos_linux64/bin/fastq-dump $file
done
#TianR. Aug 14, 2012
#Revision HMP
#for file in SRR016992.fastq SRR016993.fastq SRR016994.fastq SRR016995.fastq SRR016996.fastq SRR016997.fastq SRR016998.fastq SRR016999.fastq
###!!!!Disk Quota curse!!!!!
for file in `ls *tq`
do
        echo $file
        ###!!index is hg19####
        bowtie -a --best --strata -S -v 2 -m 3 /mnt/Storage/data/Bowtie/hg19 $file $file".SAM"
        python /mnt/Storage/home/tianr/Ich_nur/Data_Raw/ScriptsRunLog/samtobed.py $file".SAM"
        #rm $file".SAM"
	gzip $file
done
~      
