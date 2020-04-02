#!/bin/bash

#mini script to run ms for constructing subsequent sfs based on dadi parameters!
#p0=nu1, nu2,  hrf, m12, m21, me12, me21, Ts, Tsc, P, Q = params                        
#here replace the value by the estimated parameters!
model=$1

#### DEFINE DADI PARAMETER VALUES #############
nu1=3.79081971
nu2=5.27182934
b1=0.10047173
b2=0.21163122
hrf=0.0885521
m12=43.33711831
m21=20.63769348
me12=1.01587538
me21=0.19402431
Ts=2.09085953
Tsc=0.16105228
P=0.65174253
Q=0.49943849
#=-11164.004164128995
theta=1883.7033051085953
#=22354.00832825799

###CONVERRT SOME VALUE INTO MS FORMAT ####################"
m12=$(echo $(bc -l <<< " ${m12} * 2 " ))
m21=$(echo $(bc -l <<<" ${m21} * 2 "))
me12=$(echo $(bc -l <<<" ${me12} * 2 "))
me21=$(echo $(bc -l <<<" ${me21} * 2 "))
Ts=$(echo $(bc -l <<<" ${Ts} / 2 "))
Tsc=$(echo $(bc -l <<<" ${Tsc} / 2 "))
nub1=$(echo $(bc -l <<<" $nu1 * $b1 " ))
nub2=$(echo $(bc -l <<<" $nu2 * $b2 " ))

#### THEN RUN MS ACCORDING TO BEST MODEL:
if [ $model = "SC2N2m" ] ;
then 
ms 160 1 -t $theta -r 0.003 5241982 -I 2 80 80 0 \
        -m 1 2 $m12  \
        -m 2 1 $m21 \
        -n 1 $nu1 \
        -n 2 $nu2 \
        -eM $Tsc 0 \
        -ej $Ts 2 1
else
#[[ $model= "SC2N2mG" ]]  
ms 160 1 -t $theta -r 0.003 5241982 -I 2 80 80 0 \
        -m 1 2 $m12  \
        -m 2 1 $m21 \
        -n 1 $nub1 \
        -n 2 $nub2 \
	-en $nu1 1 $Tsc \
	-en $nu2 2 $Tsc \
        -eM $Tsc 0 \
        -ej $Ts 2 1
fi

exit
#fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI                        
