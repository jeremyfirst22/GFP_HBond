#!/bin/bash


molecList="
GFP_WT
GFP_F114X
GFP_F145X
GFP_F165X
GFP_M218X
GFP_N212X
GFP_Y143X
GFP_Y149X
GFP_Y92X 
GFP_WT"

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
    if [ ! -f rmsd/back_rmsd.xvg ] ; then 
        cd rmsd 
        echo 'r CRO & ! a H*' > selection.dat 
        echo 'q' >> selection.dat 
        cat selection.dat | gmx make_ndx -f ../Production/$molec.production.tpr -o cro.ndx 
        check cro.ndx 
        
        echo '18 18' | gmx rms -s ../Production/$molec.production.tpr -f ../Production/$molec.production.nopbc.xtc -n cro.ndx -o cro_alone.xvg
        check cro_alone.xvg 

        echo '4 18' | gmx rms -s ../Production/$molec.production.tpr -f ../Production/$molec.production.nopbc.xtc -n cro.ndx -o cro_rmsd.xvg
        check cro_rmsd.xvg 

        echo '4 4' | gmx rms -s ../Production/$molec.production.tpr -f ../Production/$molec.production.nopbc.xtc -o back_rmsd.xvg
        check back_rmsd.xvg 
        cd ../
    fi 
} 

r_gyrate(){
    if [ ! -d gyrate ] ; then mkdir gyrate ; fi 
    if [ ! -f gyrate/gyrate.xvg ] ; then 
        cd gyrate 
        echo '4' | gmx gyrate -s ../Production/$molec.production.tpr -f ../Production/$molec.production.xtc -o gyrate.xvg
        check gyrate.xvg
        cd ../
        fi 
} 

cro_hist147_distance(){
    if [ ! -d cro_hist147 ] ; then mkdir cro_hist147 ; fi 
    if [ ! -f cro_hist147/cro_hist147.xvg ] ; then 
        cd cro_hist147 
        echo 'r HIS & ri 147' > selection.dat 
        echo 'r CRO' >> selection.dat 
        echo 'q' >> selection.dat 
        
        hisCG=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep 147 | grep CG | awk '{print $3}'`
        hisND1=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep 147 | grep ND1 | awk '{print $3}'`
        hisCE1=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep 147 | grep CE1 | awk '{print $3}'`
        hisNE2=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep 147 | grep NE2 | awk '{print $3}'`
        hisCD2=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep 147 | grep CD2 | awk '{print $3}'`

        echo "[ his147 ] " > index.ndx
        echo "$hisCG $hisND1 $hisCE1 $hisNE2 $hisCD2" >> index.ndx 

        croCG2=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CG2 | awk '{print $3}'` 
        croCD1=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CD1 | awk '{print $3}'` 
        croCD2=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CD2 | awk '{print $3}'` 
        croCE1=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CE1 | awk '{print $3}'` 
        croCE2=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CE2 | awk '{print $3}'` 
        croCZ=`cat ../Production/$molec.npt_relax.gro | grep 66CRO | grep 66 | grep CZ | awk '{print $3}'` 

        echo "[ cro66 ] " >> index.ndx 
        echo "$croCG2 $croCD1 $croCD2 $croCE1 $croCE2 $croCZ" >> index.ndx 
         

        gmx distance -n index.ndx -s ../Production/$molec.production.tpr -f ../Production/$molec.production.xtc -select 'com of group "his147" plus com of group "cro66"' -oav cro_hist147.xvg 
        check cro_hist147.xvg
        cd ../
        fi 
} 

his147_chi1(){
    if [ ! -d his147_chi1 ] ; then mkdir his147_chi1 ; fi 
    if [ ! -f his147_chi1/his147_chi1.xvg ] ; then 
        cd his147_chi1

        hisN=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep " N " | awk '{print $3}'`
        hisCA=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep " CA " | awk '{print $3}'`
        hisCB=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep " CB " | awk '{print $3}'`
        hisCG=`cat ../Production/$molec.npt_relax.gro | grep 147HIS | grep " CG " | awk '{print $3}'`
        
        echo "[ his147_chi1 ] " > index.ndx
        echo "$hisN $hisCA $hisCB $hisCG" >> index.ndx 
        
        gmx angle -f ../Production/$molec.production.xtc -type dihedral -n index.ndx -ov his147_chi1.xvg -od his147_dist.xvg 
        check his147_chi1.xvg
        cd ../
        fi 
}

nit_chr_dist(){
    if [ ! -d nit_chr ] ; then mkdir nit_chr ; fi 
    if [ ! -f nit_chr/nit_chr.xvg ] ; then 
        cd nit_chr
        
        echo 'r CRO or CNF' > selection.dat 
        echo 'q' >> selection.dat 
        cat selection.dat | gmx make_ndx -f ../Production/$molec.production.tpr -o cnf_cro.ndx
        check cnf_cro.ndx

        echo '18' | gmx distance -s ../Production/$molec.production.tpr -f ../Production/$molec.production.xtc -n cnf_cro.ndx -oav nit_chr.xvg 
        check nit_chr.xvg
        cd ../
        fi 
}



for molec in $molecList ; do 
    printf "$molec\n\n" 
    
    cd $molec

    rmsd
    r_gyrate
    #b_b_distance
    cro_hist147_distance 
    his147_chi1
    if [ "$molec" != "GFP_WT" ] ; then 
        nit_chr_dist
        fi 
    cd ../

    done 
