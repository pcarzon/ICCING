#!/bin/bash

echo >  "$DataDir"/"$DataFile".dat

for ((m="$Begin"; m<="$End"; m++)) 
do 
   sed -i '$d' "$DataDir"/"$m"/"$DataFile".dat
   cat "$DataDir"/"$m"/"$DataFile".dat >>  "$DataDir"/"$DataFile".dat
done

 tail -n +2 "$DataDir"/"$DataFile".dat > "$DataDir"/"$DataFile".tmp && mv "$DataDir"/"$DataFile".tmp "$DataDir"/"$DataFile".dat
