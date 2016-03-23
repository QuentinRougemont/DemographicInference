#!/bin/bash

for i in `seq 1 1 20`; do

	python script_inference_anneal2_newton_mis.py -f /home/ubuntu/Script/dadi-1.6.3_modif/Lampetra&Demo/DataFolded -y LF -x LP -p 30,40,50 -f AA_28_28_indep.fs -m SI,IM,AM,SC,SC2M,AM2M, IM2M -z -l -v
	
((i++))
done


