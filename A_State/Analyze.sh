#!/bin/bash


molecList="
GFP_F114X
GFP_F145X
GFP_F165X
GFP_M218X
GFP_N212X
GFP_Y143X
GFP_Y149X
GFP_Y92X 
"

check(){
    for var in $@ ; do 
        if [ ! -s $var ] ; then 
            echo "ERROR: $var missing. " 
            exit 
            fi 
        done 
} 

rmsd(){
    if [ ! -d rmsd ] ; then mkdir rmsd ; fi 
    if [ ! -f ${molec}_rms.xvg ] ; then 
        cd rmsd 
        if [ ! -f crystal.ndx ] ; then 
            echo '4 && ri 3-230' > selection.dat  
            echo 'q' >> selection.dat 

            cat selection.dat | gmx make_ndx -f ../Production/$molec.production.tpr -o crystal.ndx 
            check crystal.ndx 
        fi 
        echo '22 22' | gmx rms -s ../Production/$molec.production.tpr -f ../Production/$molec.production.nopbc.xtc -o ${molec}_rms.xvg -n crystal.ndx 
        check ${molec}_rms.xvg 
        cd ../
    fi 
} 

r_gyrate(){
    if [ ! -d gyrate ] ; then mkdir gyrate ; fi 
    if [ ! -f gyrate/${molec}_gyrate.xvg ] ; then 
        cd gyrate 
        if [ ! -f cyrstal.ndx ] ; then 
            echo '4 && ri 3-230' > selection.dat 
            echo 'q' >> selection.dat 
            
            cat selection.dat | gmx make_ndx -f ../Production/$molec.production.tpr -o crystal.ndx 
            check crystal.ndx  
        fi 

        echo '22' | gmx gyrate -s ../Production/$molec.production.tpr -f ../Production/$molec.production.xtc -o ${molec}_gyrate.xvg -n crystal.ndx 
        check ${molec}_gyrate.xvg
        cd ../
    fi 
} 

for molec in $molecList ; do 
    printf "$molec\n\n" 
    
    cd $molec

    rmsd
    r_gyrate
    cd ../

    done 
