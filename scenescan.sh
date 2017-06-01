#!/bin/bash


# Define runname

echo 'What do you want to call your run?'

read runname

cd ~/SCENE


if [ ! -d "$runname" ]

then mkdir "Data"/"$runname"
fi



# Boundary file

echo "which device are you looking at"

read device

bdy='bdy'
dat='dat'

if [  -d boundaries/"$device" ] 

then 

cp boundaries/"$device"/"$device.$bdy" "Data"/"$runname"/bdy.txt
cp boundaries/"$device"/"$device.$dat" "Data"/"$runname"/"$device.$dat"

else

    echo 'No such device'
    exit

fi



# Parameter to be scanned

#echo "Which paramater do you want to change?"


#read param


#echo 'What values do you want to have for it'

#read -a values


echo 'Filename with parameter scans'

read filename

file="$filename.$dat"

echo $file


if [ -f "$file" ]
then
    read -r param<"$file"

    values=$(sed -n "2 p" "$file")
    
    echo "Param is $param"
    echo "Values are $values"
fi
  


cd "Data"/"$runname"

#for (( i=1; i<=${#values[@]}; i++ ))
for i in ${values[@]}
  do

      dir=$param'=_'$i
      echo "$dir"

  mkdir "$dir"

  cp bdy.txt "$dir"
  cp "$device.$dat" "$dir"

  cd "$dir"

  sed -i -e "s/'$param.*/'$param' $i/g" $device'.'$dat


  yes $device | ~/SCENEv2/trunk/scene

  yes $device | ipython ~/SCENEv2/graphs/graphs.py
 
  cd ..

done




