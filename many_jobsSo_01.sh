#!/bin/bash

###NUM=$(cat goodroot_data_lhc13cpass2_full.list | wc -l)
#####NUM=$(cat goodroot_data_lhc15n_pass1_pp502tev_7runs.list | wc -l)
#echo $NUM

for ((INDEX=0; INDEX<1340; INDEX++)) ## 1340 for EPOSLHC 355  for pythiaMonash, 1122 for perugia long 1107 per0short

do
#Crear directorio donde dejar los logs de PBS
   mkdir /home/hectorbm/testES/data1/$INDEX
   qsub -v JOBID=$INDEX anaSpherocity_01.pbs
    STATE=$?
   sleep 1
done
