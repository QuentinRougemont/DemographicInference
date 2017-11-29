#!/bin/bash


#echo "#!/bin/bash" > dadi.qarray
#echo "#$ -S /bin/bash" >> dadi.qarray
#echo "#$ -m bea" >> dadi.qarray
#echo "#$ -V" >> dadi.qarray
#echo "#$ -cwd" >> dadi.qarray
#echo "#$ -M quentin.rougemont@rennes.inra.fr" >> dadi.qarray
#echo "#$ -e /home/qrougemont/work/RAD/DADI/log/" >> dadi.qarray

echo "module load bioinfo/dadi-1.6.3_modif" > dadi.qarray


for file in *.fs ;do for i in SI AM IM SC SC2M AM2M IM2M; do for j in `seq 1 25`; do  echo "python /usr/local/bioinfo/src/dadi/dadi_modif/dadi-1.6.3_modif/inference/script_inference_anneal2_newton_mis_new_models_folded.py \
	-o $file.$i.$j  -y LF -x LP -p 30,40,50 -f $file -m $i -l -v -z &>> $file_$i.log " >> dadi.qarray ; done ; done ; done

#script_inference_anneal2_newton_mis_new_models_folded.py 


#for file in SPECTRUM2/* ; do for i in  AM IM SC SC2M AM2M IM2M ; do qsub -q workq -v INPUT=$file  dadi.qarray ; done; done

