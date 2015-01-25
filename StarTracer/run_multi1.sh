#!/bin/bash

start=$1
end=$2
grpmask="ST"
snapmask="snap"
DirFromSnap="/ptmp1/h009z/h009zac/BOX64/SNAPS/"

DirFromGrp="/ptmp1/h009z/h009zac/BOX64/FOF1.0/"

wdir="/ptmp1/h009z/h009zac/BOX64/FOF1.0/"
fendgrp="_.grp.DAT"


while [ $start -ge "$end" ]
do
  if [ "$start" -le "9" ]
      then
      snap=00$start
  else
      if [ "$start" -le "99" ]
          then
          snap=0$start
      else
          if [ "$start" -le "999" ]
              then
              snap=$start
          fi
      fi
  fi
start2=$(($start-2))

  if [ "$start2" -le "9" ]
      then
      snap2=00$start2
  else
      if [ "$start2" -le "99" ]
          then
          snap2=0$start2
      else
          if [ "$start2" -le "999" ]
              then
              snap2=$start2
          fi
      fi
  fi
  
  fsnap1=$DirFromSnap$snapmask"_"$snap  
  fsnap2=$DirFromSnap$snapmask"_"$snap2  
  cat1=$DirFromGrp$grpmask"_"$snap$fendgrp 
  cat2=$DirFromGrp$grpmask"_"$snap2$fendgrp 
  echo $fsnap1
  echo $fsnap2
  echo $cat1
  echo $cat2
  
  if [ -e $fsnap1  ] && [ -e $fsnap2  ]&& [ -e $cat1  ] && [ -e $fcat2  ]
      then
########################  
      ~/bin/StarTracer.x $fsnap1 $fsnap1 $cat1 $cat2     
      #./mysql -uroot < $DirFrom"ST_"${snap}"_.grp.DAT.mysql"   
  fi
########################
start=$(($start-1))
#exit
done


