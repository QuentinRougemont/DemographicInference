#!/bin/bash

#mini script to run ms for constructing subsequent sfs based on dadi parameters!
#p0=nu1, nu2,  hrf, m12, m21, me12, me21, Ts, Tsc, P, Q = params                        
#p0=nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, P, Q = params
#here replace the value by the estimated parameters!
model=$1

#### DEFINE DADI PARAMETER VALUES #############
nu1=0.142100567
nu2=99.0451665
b1=12.2705553
b2=0.0262736942
hrf=0.100098352
m12=13.5271892
m21=0.886416416
me12=1.34613939
me21=0.0182086547
Ts=1.48643243
P=0.650403599
Q=0.205321952
#=-9339.78909169245
theta=913.1069130897766
AIC=18703.5781833849

###CONVERRT SOME VALUE INTO MS FORMAT ####################"
m12=$(echo $(bc -l <<< " ${m12} * 2 " ))
m21=$(echo $(bc -l <<<" ${m21} * 2 "))
me12=$(echo $(bc -l <<<" ${me12} * 2 "))
me21=$(echo $(bc -l <<<" ${me21} * 2 "))
Ts=$(echo $(bc -l <<<" ${Ts} / 2 "))
#Tsc=$(echo $(bc -l <<<" ${Tsc} / 2 "))
nub1=$(echo $(bc -l <<<" $nu1 * $b1 " ))
nub2=$(echo $(bc -l <<<" $nu2 * $b2 " ))
nu1hrf=$(echo $(bc -l <<<"$nu1 * hrf"))
nu2hrf=$(echo $(bc -l <<<"$nu2 * hrf"))
bnu1hrf=$(echo $(bc -l <<<"$nub1 * hrf"))
bnu2hrf=$(echo $(bc -l <<<"$nub2 * hrf"))

#### THEN RUN MS ACCORDING TO BEST MODEL:
if [ $model = "IM2N2m" ] ;
then 
ms 160 1 -t $theta -r 0.003 1965875 -I 2 80 80 0 \
        -m 1 2 $m12  \
        -m 2 1 $m21 \
        -n 1 $nu1hrf \
        -n 2 $nu2hrf \
        -ej $Ts 2 1
else

#[[ $model= "IM2N2mG" ]]  
ms 160 1 -t $theta -r 0.003 1965875 -I 2 80 80 0 \
        -m 1 2 $m12  \
        -m 2 1 $m21 \
        -n 1 $bnu1hrf\
        -n 2 $bnu2hrf\
	-en $nu1hrf 1 $Ts \
	-en $nu2hrf 2 $Ts \
        -ej $Ts 2 1
fi

exit
#fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI                        
