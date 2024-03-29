#!/bin/bash

for file in *.cif
do
    mkdir ${file%%.*}
    cp *.xyz gas.xyz 
    cp $file gas.xyz ${file%%.*}
    cd ${file%%.*}
    cp ~/Documents/Arbalign/bin/* ./
    rm run_arba.sh
    sed -i "s/CIFFF/$file/g" run_arbalign.py
    sed -i "s/NNN/gas.xyz/g" run_arbalign.py
    sed -i "s/FINN/final/g" run_arbalign.py
    echo "What is Z for this crystal? :"
    read Zeta
    sed -i "s/PPPPP/$Zeta/g" run_arbalign.py
    echo "Which axis [A,B,C] does this crystal dimerize on?: "
    read axis
    if [ $axis == "a" ] || [ $axis == "A" ];
        then
        sed -i "s/xxx/2/g" run_arbalign.py 
        sed -i "s/yyy/1/g" run_arbalign.py
        sed -i "s/zzz/1/g" run_arbalign.py
    elif [ $axis == "b" ] || [ $axis == "B" ];
        then
        sed -i "s/xxx/1/g" run_arbalign.py 
        sed -i "s/yyy/2/g" run_arbalign.py
        sed -i "s/zzz/1/g" run_arbalign.py
    elif [ $axis == "c" ] || [ $axis == "C" ];
        then
        sed -i "s/xxx/1/g" run_arbalign.py 
        sed -i "s/yyy/1/g" run_arbalign.py
        sed -i "s/zzz/2/g" run_arbalign.py
    else
        echo "That is not a suitable answer"
    fi
python3 run_arbalign.py
done   