#!/bin/bash


# Define runname

echo 'What do you want to call your run?'

read runname

cd ~/SCENEv2/


if [ ! -d "Data"/"$runname" ]

then mkdir "Data"/"$runname"
fi



# Boundary file

echo "which device are you looking at"

read device

bdy='bdy'
dat='dat'

if [  -d boundaries/"$device" ] 

then 

cp boundaries/"$device"/bdy.txt "Data"/"$runname"/bdy.txt
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

    param1=$(sed -n "1p" "$file")    
    val1=$(sed -n "2 p" "$file")
    param2=$(sed -n "3p" "$file")
    val2=$(sed -n "4 p" "$file")
   
    echo "Param 1 is $param1"
    echo "Values 1 are $val1"
    echo "param 2 is $param2"
    echo "Values 2 are $val2"
else
    echo 'Parameter data not found'
    exit

fi
  


cd "Data"/"$runname"
count=0

#for (( i=1; i<=${#values[@]}; i++ ))
for i in ${val1[@]}
do

    sed -i -e "s/'$param1.*/'$param1' $i/g" $device'.'$dat
    for j in ${val2[@]}
	     
    do

	count=$(($count+1))
	echo "$count"
	dir=$param1'='$i$param2'='$j
	echo "$dir"

	mkdir "$dir"

	#Copy data from head directory
	cp bdy.txt "$dir"

	cp "$device.$dat" "$dir"

	if [ $count != 1 ]
	then
        
	    cp fort.7 "$dir"

	fi


	
	cd "$dir"


	sed -i -e "s/'$param2.*/'$param2 ' $j/g" $device'.'$dat

	echo "Run for $param1 = $i and $param2 = $j" 
	yes $device | ~/SCENEv2/trunk/scene

	yes $device | ipython ~/SCENEv2/graphs/graphs.py

	#Copy restart data to above directory
	cp fort.7 ../
	
	cd ..

	#If its the first run, change input file so next run start from prev eqbm
	if [ $count == 1 ]
	then
	    echo 'Setting input to restart'
	    sed -i -e "s/'icon'*/'icon' -1/g" $device'.'$dat
	fi
	
    done
    
done
