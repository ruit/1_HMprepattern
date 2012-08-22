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
