#!/bin/bash
#Tian R. & Zhang Y., Dec.8, 2011
#Jan 05, 2012
# For single HMs
# For single HM

for HM in H2AZ_Pro H3K27me1_Body H3K27me3_Pro H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Pro H4K20me1_Body 
	do
		cat 2012.01.5_singleHM.R | sed "s/H2AZ_Pro/$HM/g" > temp.R
		Rscript temp.R
	done

#rm temp.R

#########################
#For double HMs
#for HMa in H2AZ_Body H3K27me1_Body H3K27me3_Body H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Body H4K20me1_Pro
#	do
#		for HMb in H2AZ_Body H3K27me1_Body H3K27me3_Body H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Body H4K20me1_Pro 
#			do
#				if [ $HMa == $HMb ]; then
#					echo "the same HMs: $HMa!"
#				else
#					cat 12.8_subset.features.R | sed "s/H2AZ_Body1/$HMa/g" | sed "s/H2AZ_Body2/$HMb/g" > Temp.R
#					Rscript Temp.R					
#				fi
#			done
#	done

##########################
#For tripple HMs
#for HMA in H2AZ_Body H3K27me1_Body H3K27me3_Body H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Body H4K20me1_Pro
 #      do
  #             for HMB in H2AZ_Body H3K27me1_Body H3K27me3_Body H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Body H4K20me1_Pro 
   #                    do
    #                           if [ $HMA == $HMB ]; then
     #                                  echo "the same HMs: $HMA!"
      #                         else
		#			for HMC in H2AZ_Body H3K27me1_Body H3K27me3_Body H3K36me3_Body H3K4me1_Body H3K4me3_Pro H3K9me1_Body H3K9me3_Body H4K20me1_Pro 
		#				do
		#					
		#					if [ $HMA == $HMC ] || [ $HMB == $HMC ];then
		#						echo "The Same HMs: $HMC!"
		#					else
		#						cat 12.8_subset.features.R | sed "s/H2AZ_BodyA/$HMA/g" | sed "s/H2AZ_BodyB/$HMB/g" | sed "s/H2AZ_BodyC/$HMC/g" > TEMP.R 
		#						Rscript TEMP.R
		#					fi
		#				done
		#		fi
		#	done
#	done

