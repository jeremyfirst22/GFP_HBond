#!/bin/bash

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} > "
    exit 
}


if [ -z $1 ] ; then 
    usage 
    fi 

fileName=$1 
if [ ! -f $fileName ] ; then 
    echo "ERROR: $fileName not found " 
    exit 
    fi 
if [[ $fileName == *.pdb ]] ; then 
    MOLEC=$(basename $fileName)
    MOLEC=${MOLEC%.*}
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -d $MOLEC ] ; then mkdir $MOLEC ; fi 
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/ ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log 
errFile=$TOP/$MOLEC/$MOLEC.err 
FF=$TOP/GMXFF

check(){
   for var in $@ ; do 
        if [ ! -s $var ] ; then 
            echo ; echo $var missing, exitting... 
            exit 
            fi 
        done 
}

clean(){
    if [ -d amber03.ff ] ; then rm -r amber03.ff *.dat ; fi 
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi 

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi 
    
    if [ ! -d $dirName/amber03.ff ] ; then 
        if [ -d $FF/amber03.ff ] ; then 
            cp -r $FF/amber03.ff $dirName/amber03.ff 
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found. Standard ff not yet combined with CRO augmented force field." 
            exit 
            fi 
        fi 
}

protein_min(){
    printf "\t\tProtein STEEP................." 
    if [ ! -f Protein_min/$MOLEC.solute_min.gro ] ; then 
        create_dir Protein_min

        cp $MOLEC.pdb Protein_min
        cd Protein_min
        
        gmx pdb2gmx -f $MOLEC.pdb -o $MOLEC.gro -p $MOLEC.top -ff amber03 -water tip3p >> $logFile 2>> $errFile
        check $MOLEC.top $MOLEC.gro 
        
        gmx grompp -f $MDP/solute_steep.mdp -c $MOLEC.gro -p $MOLEC.top -o $MOLEC.solute_min.tpr >> $logFile 2>> $errFile
        check $MOLEC.solute_min.tpr 

        gmx mdrun -deffnm $MOLEC.solute_min >> $logFile 2>> $errFile
        check $MOLEC.solute_min.gro 
        
        gmx editconf -f $MOLEC.solute_min.gro -o $MOLEC.min.pdb >> $logFile 2>> $errFile
        check $MOLEC.min.pdb 
        cp $MOLEC.min.pdb ../

        clean
        printf "Success\n" 
        cd ../

    else 
        printf "Skipped\n" 
        fi 
    check Protein_min/$MOLEC.solute_min.gro 
}

solvate(){
    printf "\t\tSolvating....................." 
    if [ ! -f Solvate/$MOLEC.neutral.gro ] ; then 
        create_dir Solvate

        cp Protein_min/$MOLEC.solute_min.gro Solvate
        cd Solvate

        gmx pdb2gmx -f $MOLEC.solute_min.gro -o $MOLEC.gro -p $MOLEC.top -water tip3p -ff amber03 >> $logFile 2>> $errFile 
        check $MOLEC.gro $MOLEC.top 
        
        gmx editconf -f $MOLEC.gro -bt octahedron -box 8 -o temp.centered.gro >> $logFile 2>> $errFile 
        check temp.centered.gro 
        
        gmx solvate -cp temp.centered.gro -o temp.solvated.gro >> $logFile 2>> $errFile  
        check temp.solvated.gro 
        
        gmx pdb2gmx -f temp.solvated.gro -water tip3p -ff amber03 -o temp.$MOLEC.solvated.gro -p temp.$MOLEC.solvated.top >> $logFile 2>> $errFile 
        check temp.$MOLEC.solvated.gro temp.$MOLEC.solvated.top

        gmx grompp -f $MDP/vac_md.mdp -p temp.$MOLEC.solvated.top -c temp.$MOLEC.solvated.gro -o temp.genion.tpr >> $logFile 2>> $errFile 
        check temp.genion.tpr

        echo SOL | gmx genion -s temp.genion.tpr -neutral -nname 'Cl-' -pname 'Na+' -o temp.neutral.gro >> $logFile 2>> $errFile 
        check temp.neutral.gro 

        gmx pdb2gmx -f temp.neutral.gro -water tip3p -ff amber03 -p $MOLEC.neutral.top -o $MOLEC.neutral.gro >> $logFile 2>> $errFile 
        check $MOLEC.neutral.gro 

        clean
        printf "Success\n" 
        cd ../

   else 
       printf "Skipped\n"      
       fi 
   check Solvate/$MOLEC.neutral.gro Solvate/$MOLEC.neutral.top 
}

solvent_min(){
    printf "\t\tSolvent minimization.........." 
    if [ ! -f Solvent_min/$MOLEC.npt_relax.gro ] ; then 
        if [ ! -f Solvent_min/$MOLEC.minimize.gro ] ; then 
            create_dir Solvent_min

            cp Solvate/$MOLEC.neutral.gro Solvent_min
            cp Solvate/$MOLEC.neutral.top Solvent_min
            cp Solvate/*.itp Solvent_min
            cd Solvent_min

            gmx grompp -f $MDP/solvent_min.mdp -c $MOLEC.neutral.gro -p $MOLEC.neutral.top -o $MOLEC.minimize.tpr >> $logFile 2>> $errFile 
            check $MOLEC.minimize.tpr

            gmx mdrun -deffnm $MOLEC.minimize >> $logFile 2>> $errFile 
            check $MOLEC.minimize.gro
            fi 
            
        if [ ! -f $MOLEC.nvt_relax.gro ] ; then 
            gmx grompp -f $MDP/solvent_nvt_relax.mdp -c $MOLEC.minimize.gro -p $MOLEC.neutral.top -o $MOLEC.nvt_relax.tpr >> $logFile 2>> $errFile 
            check $MOLEC.nvt_relax.tpr
            
            gmx mdrun -deffnm $MOLEC.nvt_relax >> $logFile 2>> $errFile 
            check $MOLEC.nvt_relax.gro 
            fi 

        if [ ! -f $MOLEC.npt_relax.gro ] ; then 
            gmx grompp -f $MDP/solvent_npt_relax.mdp -c $MOLEC.nvt_relax.gro -p $MOLEC.neutral.top -o $MOLEC.npt_relax.tpr >> $logFile 2>> $errFile 
            check $MOLEC.npt_relax.tpr 

            gmx mdrun -deffnm $MOLEC.npt_relax >> $logFile 2>> $errFile 
            check $MOLEC.npt_relax.gro 
            fi 

        clean
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
    check Solvent_min/$MOLEC.npt_relax.gro Solvent_min/$MOLEC.neutral.top 
} 

production_run(){
    printf "\t\tProduction run ..............." 
    if [ ! -f Production/$MOLEC.production.nopbc.gro ] ; then 
        create_dir Production
        cp Solvent_min/$MOLEC.npt_relax.gro Production
        cp Solvent_min/$MOLEC.neutral.top Production
        cp Solvent_min/*.itp Production
        cd Production
        
        if [ ! -f $MOLEC.production.gro ] ; then 
            if [ ! -f $MOLEC.production.tpr ] ; then 
                gmx grompp -f $MDP/production_gfp.mdp -o $MOLEC.production.tpr -p $MOLEC.neutral.top -c $MOLEC.npt_relax.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.production.tpr 
            if [ -f $MOLEC.production.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC.production -cpi $MOLEC.production.cpt >> $logFile 2>> $errFile 
            else 
                gmx mdrun -deffnm $MOLEC.production >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.production.gro 
        
        if [ ! -f $MOLEC.production.nopbc.gro ] ; then 
            echo '1 0' | gmx trjconv -f $MOLEC.production.xtc -s $MOLEC.production.tpr -pbc mol -ur compact -center -o $MOLEC.production.nopbc.xtc >> $logFile 2>> $errFile 
            echo '1 0' | gmx trjconv -f $MOLEC.production.gro -s $MOLEC.production.tpr -pbc mol -ur compact -center -o $MOLEC.production.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.production.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
}

analyze_hbond(){
    printf "\t\tAnalyzing Hbond content ......" 
    if [[ ! -f HBond/cnf_num.xvg || ! -f HBond/cro_num.xvg ]] ; then 
        create_dir HBond
        cd HBond
    
        if [ ! -f cnf_num.xvg ] ; then 
            echo "r CNF & a NH" > selection.dat 
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.gro -o cnf_nh.ndx >> $logFile 2>> $errFile 
            check cnf_nh.ndx 

            echo '18 14 18' | gmx hbond -f ../Production/$MOLEC.production.xtc -s ../Production/$MOLEC.production.tpr -n cnf_nh.ndx -shell 1 -r 0.3 -a 20 -num cnf_num.xvg >> $logFile 2>> $errFile 
            fi 
        check cnf_num.xvg 
            
        if [ ! -f cro_num.xvg ] ; then 
            echo "r CROn & a OH or a HO" > selection.dat 
            echo "r CROn & a OH" >> selection.dat
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.gro -o cro_o_h.ndx >> $logFile 2>> $errFile
            check cro_o_h.ndx

            echo '18 14 19' | gmx hbond -f ../Production/$MOLEC.production.xtc -s ../Production/$MOLEC.production.tpr -n cro_o_h.ndx -shell 1 -r 0.3 -a 20 -num cro_num.xvg >> $logFile 2>> $errFile 
            fi 
        check cro_num.xvg

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

analyze_hbond_nit(){
    printf "\t\tAnaylzing HBond_nit content..."
    if [ ! -s HBond_nit/$MOLEC.hb_count.xvg ] ; then 
        create_dir HBond_nit
        cd HBond_nit

    if [ ! -d ../Production/amber03.ff ] ; then 
        cp $FF/*.dat ../Production/. 
        cp -r $FF/amber03.ff ../Production/.
        fi 
    check ../Production/amber03.ff/forcefield.itp 

        ## We use veriosn 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
        ## We allow for warnings, since we are generated .tpr from a gromacs 5 mdp file. We are only inserting
        ## atoms this should not matter. 
        if [ ! -f $MOLEC.production.v4.tpr ] ; then 
            grompp -f $MDP/production_gfp.mdp -o $MOLEC.production.v4.tpr -p ../Production/$MOLEC.neutral.top -c ../Production/$MOLEC.npt_relax.gro -maxwarn 3 >> $logFile 2>> $errFile 
            fi
        check $MOLEC.production.v4.tpr 
    
        CT=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep NH | awk '{print $3}'`
        #echo $CT $NH
        
        if [ ! -s $MOLEC.hb_count.xvg ] ; then  
            $HOME/andrews_gmx/g_nitrile_hbond/g_nitrile_hbond \
                -s $MOLEC.production.v4.tpr \
                -f ../Production/$MOLEC.production.nopbc.xtc \
                -a1 $CT -a2 $NH \
                -select "resname SOL and same residue as within 0.5 of resname CNF and name NH" \
                -o $MOLEC.frame_hb.xvg -op $MOLEC.persistent.xvg \
                -oa $MOLEC.hb_count.xvg -or $MOLEC.geometry.xvg  >> $logFile 2>> $errFile
            check $MOLEC.hb_count.xvg $MOLEC.geometry.xvg $MOLEC.persistent.xvg $MOLEC.frame_hb.xvg 
        fi 
        
        check $MOLEC.hb_count.xvg 
        printf "Success\n" 
        clean
        cd ../
    else 
        printf "Skipped\n" 
        fi 
    check HBond_nit/$MOLEC.hb_count.xvg 
}

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
protein_min
solvate
solvent_min
production_run 
analyze_hbond
analyze_hbond_nit
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
