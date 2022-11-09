#!/bin/bash

for file in *.cif
do
    mkdir ${file%%.*}
    mv *.xyz gas.xyz 
    mv $file gas.xyz ${file%%.*}
    cd ${file%%.*}
    cp ~/Documents/Arbalign/bin/* ./
    sed -i "s/CIFFF/$file/g" run_arbalign.py
    sed -i "s/NNN/gas.xyz/g" run_arbalign.py
    sed -i "s/FINNN/${file%%.*}/g" run_arbalign.py
    echo "Which axis [A,B,C] does this crystal dimerize on?: "
    read axis
    if [$axis == "a"] || [$axis == "A"]
        then
        sed -i "s/xxx/2/g" run_arbalign.py 
        sed -i "s/yyy/1/g" run_arbalign.py
        sed -i "s/zzz/1/g" run_arbalign.py
    elif [$axis == "b"] || [$axis == "B"]
        then
        sed -i "s/xxx/1/g" run_arbalign.py 
        sed -i "s/yyy/2/g" run_arbalign.py
        sed -i "s/zzz/1/g" run_arbalign.py
    elif [$axis == "b"] || [$axis == "B"]
        then
        sed -i "s/xxx/1/g" run_arbalign.py 
        sed -i "s/yyy/1/g" run_arbalign.py
        sed -i "s/zzz/2/g" run_arbalign.py
    else
        echo "That is not a suitable answer"
    fi
done   
