#! /bin/bash

for j in $(ls ./x*); do
	for i in $(cat $j); do
		
		echo "$i"
		echo "split"
		split=$(echo $i | sed -e 's/+/ /g')
		echo "$split"
		echo "waffle"
		echo "waffle"


	done
        wait
done
