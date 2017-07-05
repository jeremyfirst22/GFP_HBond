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
FORCE_TOOLS=/Users/jeremyfirst/force_calc_tools

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
    printf "\t\tProtein STEEP...................." 
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
    printf "\t\tSolvating........................" 
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
    printf "\t\tSolvent minimization............." 
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
    printf "\t\tProduction run..................." 
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
    printf "\t\tAnalyzing Hbond content.........." 
    if [ ! -f hbond/rates.txt ] ; then 
        create_dir hbond
        cd hbond
        
        touch empty.ndx 
        if [ ! -f cnf_nh.xvg ] ; then 
            echo "r CNF & a NH" > selection.dat 
            echo "! r CNF" >> selection.dat
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.gro -o cnf_nh.ndx -n empty.ndx >> $logFile 2>> $errFile
            check cnf_nh.ndx

            echo '0 1 0' | gmx hbond -f ../Production/$MOLEC.production.xtc -s ../Production/$MOLEC.production.tpr -n cnf_nh.ndx -shell 0.5 -r 0.3 -da -num -ac -ang -dist -don -dan > hbond.log 2>> $errFile 
            fi 
        check hbnum.xvg danum.xvg donor.xvg hbac.xvg hbang.xvg hbdist.xvg 
         
        grep -A 100 "Doing autocorrelation according to the theory of Luzar and Chandler." hbond.log > rates.txt  

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

analyze_hbond_nit(){
    printf "\t\tAnaylzing HBond_nit content......"
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

min_dist(){
    printf "\t\tCalculating min_dist to water...."
    if [ ! -s min_dist/cnf_water.xvg ] ; then 
        create_dir min_dist
        cd min_dist
    
        if [ ! -f cnf_water.ndx ] ; then 
            echo "resname CNF and name NH" > selection.dat 
            echo "group Water and (name HW1 or name HW2)" >> selection.dat 

            cat selection.dat | gmx select -s ../Production/$MOLEC.production.tpr -f ../Production/$MOLEC.production.nopbc.xtc -on cnf_water.ndx >> $logFile 2>> $errFile 
            fi 
        check cnf_water.ndx 

        echo '0 1' | gmx mindist -f ../Production/$MOLEC.production.nopbc.xtc -n cnf_water.ndx -xvg none -od cnf_water.xvg >> $logFile 2>> $errFile

        check cnf_water.xvg 
        printf "Success\n" 
        clean 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
} 

force_calc(){
    printf "\n\t\tCalculating force:\n" 
    if [[ ! -f force_calc/$MOLEC.solvent_rxn_field.projected.xvg || ! -f force_calc/$MOLEC.external_field.projected.xvg || ! -f force_calc/$MOLEC.total_field.projected.xvg ]] ; then 

    if [ ! -f $FORCE_TOOLS/g_insert_dummy_atom ] ; then 
        printf "\t\t\tERROR: Force tools not found. Skipping force calc\n" 
        return  
        fi 

    create_dir force_calc
    cd force_calc 

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
    
    printf "\t\t\tInserting dummy atoms............................" 
    if [ ! -s $MOLEC.with_dummy.xtc ] ; then 
        $FORCE_TOOLS/g_insert_dummy_atom -s $MOLEC.production.v4.tpr -f ../Production/$MOLEC.production.nopbc.xtc -o $MOLEC.with_dummy.xtc -a1 $CT -a2 $NH >> $logFile 2>> $errFile 
        printf "Done\n" 
    else 
        printf "Skipped\n" 
        fi 
    check $MOLEC.with_dummy.xtc

    if [ ! -s $MOLEC.with_dummy.gro ] ; then 
        $FORCE_TOOLS/g_insert_dummy_atom -s $MOLEC.production.v4.tpr -f ../Production/$MOLEC.production.nopbc.gro -o $MOLEC.with_dummy.gro -a1 $CT -a2 $NH >> $logFile 2>> $errFile 
        fi 
    check $MOLEC.with_dummy.gro 

    if [ ! -s $MOLEC.with_dummy.top ] ; then 
        if [ "${MOLEC: -1}" == "H" ] ; then 
            gmx pdb2gmx -f $MOLEC.with_dummy.gro -p $MOLEC.with_dummy.top -water tip3p -ff amber03 >> $logFile 2>> $errFile 
        else 
            gmx pdb2gmx -f $MOLEC.with_dummy.gro -p $MOLEC.with_dummy.top -water tip3p -ff amber03 >> $logFile 2>> $errFile 
            fi 
        fi 
    check $MOLEC.with_dummy.top 
    
    ##Find new atom numbers 
    CT=`grep CNF $MOLEC.with_dummy.gro | grep CT | awk '{print $3}'`
    NH=`grep CNF $MOLEC.with_dummy.gro | grep NH | awk '{print $3}'`

    echo "[ probe ]" > probe.ndx 
    echo "$CT $NH" >> probe.ndx 

    echo "[ protein ]" > protein.ndx 
    grep -v TCHG $MOLEC.with_dummy.gro | grep -v SOL | grep -v Na | grep -v Cl | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 

    cp $MOLEC.with_dummy.top $MOLEC.total_field.top 

    if [ ! -s $MOLEC.solvent_rxn_field.top ] ; then 
        $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top protein.ndx $MOLEC.solvent_rxn_field.top >> $logFile 2>> $errFile 
        fi 

    if [ ! -s $MOLEC.external_field.top ] ; then 
        $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top probe.ndx $MOLEC.external_field.top >> $logFile 2>> $errFile 
        fi 
    check $MOLEC.total_field.top $MOLEC.external_field.top $MOLEC.solvent_rxn_field.top 

    for field in total_field external_field solvent_rxn_field ; do 
        printf "\t\t%20s..." "$field" 

        ##Extract forces 
        if [ ! -s $MOLEC.$field.projected.xvg ] ; then 
            printf "forces..." 
            if [ ! -s $MOLEC.$field.xvg ] ; then 
                if [ ! -s $MOLEC.$field.tpr ] ; then 
                    gmx grompp -f $MDP/rerun.mdp -p $MOLEC.$field.top -c $MOLEC.with_dummy.gro -o $MOLEC.$field.tpr  >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.$field.tpr 
 
                if [ ! -s $MOLEC.$field.trr ] ; then 
                    gmx mdrun -rerun $MOLEC.with_dummy.xtc -s $MOLEC.$field.tpr -deffnm $MOLEC.$field >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.$field.trr 

                echo 2 | gmx traj -f $MOLEC.$field.trr -s $MOLEC.$field.tpr -of $MOLEC.$field.xvg -xvg none  >> $logFile 2>> $errFile 
                rm $MOLEC.$field.trr 
            fi 
            check $MOLEC.$field.xvg 

            ##extract postions for bond vector
            printf "positions..." 
            if [ ! -s $MOLEC.positions.xvg ] ; then 
                gmx traj -f $MOLEC.with_dummy.xtc -s $MOLEC.$field.tpr -n probe.ndx -ox $MOLEC.positions.xvg -xvg none  >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.positions.xvg 

            ##project force along bond vector 
            printf "Projecting..." 
            $FORCE_TOOLS/get_force.py $MOLEC.positions.xvg $MOLEC.$field.xvg $MOLEC.$field.projected.xvg 
            check $MOLEC.$field.projected.xvg 
            printf "Done\n"  
        else 
            printf "..................................Skipped\n" 
            fi 
        done 
    check $MOLEC.total_field.projected.xvg $MOLEC.external_field.projected.xvg $MOLEC.solvent_rxn_field.projected.xvg 
    clean 
    cd ../Production/
    clean
    cd ../
    
    else 
        printf "\t\t\t\t  ...............Skipped\n" 
        fi 
    printf "\n" 
}

force_calc_APBS(){
    tot=50000
    numFrames=100
    timeStep=$(echo "$tot / $numFrames" | bc) #ps

    printf "\t\tForce Calc APBS:\n"

    force_done=1 ##True 
    for frame in $(seq 0 $timeStep $tot) ; do 
        if [[ ! -f force_calc_APBS/time_${frame}_1.txt || ! -f force_calc_APBS/time_${frame}_78.txt ]] ; then 
            force_done=0 ##Not done 
            fi 
        done 
    if [[ ! -f force_calc_APBS/$MOLEC.coloumb_field.out || ! -f force_calc_APBS/$MOLEC.rxn_field.out ]] ; then 
            force_done=0
            fi 

    if [ ${force_done} -eq 0 ] ; then 

        if [[ $timeStep -lt 20 ]] ; then 
            echo "ERROR: Requested time step is $timeStep. Frames only printed every 20 ps"
            exit ; fi 
        if (( $tot % 4 )) ; then 
            echo "WARNING: $tot % 4 != 0. There may not be frames for requested time steps" 
            fi 

        check ../free_energy_files/AMBER.DAT ../free_energy_files/AMBER.names
        create_dir force_calc_APBS 
        cd force_calc_APBS

        cp ../../free_energy_files/AMBER.DAT . 
        cp ../../free_energy_files/AMBER.names . 
        
        printf "\t\t\tPre-Compress............."
        if [ ! -f $MOLEC.compress.xtc ] ; then 
            echo '0' | gmx trjconv -f ../Production/$MOLEC.production.nopbc.xtc -s ../Production/$MOLEC.production.tpr -o $MOLEC.compress.xtc -b 0 -e $tot -dt $timeStep -tu ps >> $logFile 2>> $errFile 
            check $MOLEC.compress.xtc 
            printf "Complete\n" 
        else 
            printf "Skipped\n" 
            fi 
        
        for frame in $(seq 0 $timeStep $tot) ; do 
            printf "\t\t\tReading %5i of %5i..." $frame $tot
            if [[ ! -f time_${frame}_78.txt || ! -f time_${frame}_1.txt ]] ; then 

                if [ ! -f time_${frame}.pdb ] ; then 
                    echo '1' | gmx trjconv -f ../Production/$MOLEC.production.nopbc.xtc -s ../Production/$MOLEC.production.tpr -o time_${frame}.pdb -dump $frame -tu ps >> $logFile 2>> $errFile 
                    check time_${frame}.pdb 
                    fi 
            
                ####For whatever reason, AOBS only support 3-character names for residues 
                sed "s/CROn/CRO /" time_${frame}.pdb > temp.pdb      
                mv temp.pdb time_${frame}.pdb 

                if ! grep -sq CRO time_${frame}.pdb ; then 
                    printf "ERROR: We left CRO behind! \n" 
                    exit
                    fi 

                if [ ! -f time_${frame}.pqr ] ; then 
                    pdb2pqr time_${frame}.pdb time_${frame}.pqr --userff AMBER.DAT --usernames AMBER.names --assign-only >> $logFile 2>> $errFile
                    fi 
                check time_${frame}.pqr 

                if [ ! -f time_${frame}+DUM.pqr ] ; then 
                    python ../../free_energy_files/add_dummy_pqr.py time_${frame}.pqr 0.1 CNF CT NH > time_${frame}+DUM.pqr 
                    fi  
                check time_${frame}+DUM.pqr 
        
                sed "s/SDIE/78/" ../../free_energy_files/force_temp.in | sed "s/FRAME/${frame}/" >> time_${frame}_78.in 
                check time_${frame}_78.in 
                #echo ${frame}

                apbs time_${frame}_78.in >> $logFile 2>> $errFile 
                check time_${frame}_78.txt 

                sed "s/SDIE/1/" ../../free_energy_files/force_temp.in | sed "s/FRAME/${frame}/" >> time_${frame}_1.in 
                check time_${frame}_1.in 

                apbs time_${frame}_1.in >> $logFile 2>> $errFile 
                check time_${frame}_1.txt

                printf "Complete\n" 
            else 
                printf "Skipped\n" 
                fi 
            done 
        
        printf "\t\t\tCompiling fields........."
        if [ -f $MOLEC.coloumb_field.out ] ; then 
            rm $MOLEC.coloumb_field.out 
            fi 
        if [ -f $MOLEC.rxn_field.out ] ; then 
            rm $MOLEC.rxn_field.out
            fi 
        for frame in $(seq 0 $timeStep $tot) ; do 
            ../../free_energy_files/read_apbs_rxn_field time_${frame}+DUM.pqr time_${frame}_78.txt time_${frame}_1.txt >> $MOLEC.rxn_field.out 2>> /dev/null 
            check $MOLEC.rxn_field.out 

            python ../../free_energy_files/analytic_coloumb.py time_${frame}.pqr CNF NH CT >> $MOLEC.coloumb_field.out 2>> /dev/null 
            check $MOLEC.coloumb_field.out
            done 
        printf "Complete\n" 


        check $MOLEC.coloumb_field.out $MOLEC.rxn_field.out 
        clean 
        cd ../
        
    else 
        printf "\t\t\t\t  ...............Skipped\n" 
        fi 
    printf "\n" 
} 

sasa(){
    printf "\t\tAnalyzing SASA..................." 
    if [[ ! -f sasa/$MOLEC.nh_cz.xvg || ! -f sasa/$MOLEC.cnf.xvg || ! -f sasa/$MOLEC.nh.xvg ]] ; then 
        create_dir sasa 
        cd sasa 
    
        if [ ! -f $MOLEC.nh_cz.xvg ] ; then 
            echo "a NH | a CZ && r CNF" > selection.dat
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.nopbc.gro -o nh_cz.ndx >> $logFile 2>> $errFile 
            check nh_cz.ndx 

            gmx sasa -s ../Production/$MOLEC.production.tpr -f ../Production/$MOLEC.production.nopbc.xtc -surface 'Protein' -output '"NH_CZ_&_CNF"' -o $MOLEC.nh_cz.xvg -n nh_cz.ndx >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nh_cz.xvg 

        if [ ! -f $MOLEC.nh.xvg ] ; then 
            echo "a NH && r CNF" > selection.dat
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.nopbc.gro -o nh.ndx >> $logFile 2>> $errFile 
            check nh.ndx 

            gmx sasa -s ../Production/$MOLEC.production.tpr -f ../Production/$MOLEC.production.nopbc.xtc -surface 'Protein' -output '"NH_&_CNF"' -o $MOLEC.nh.xvg -n nh.ndx >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nh.xvg 
            
        if [ ! -f $MOLEC.cnf.xvg ] ; then 
            gmx sasa -s ../Production/$MOLEC.production.tpr -f ../Production/$MOLEC.production.nopbc.xtc -o $MOLEC.cnf.xvg -surface 'Protein' -output 'resname CNF' >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.cnf.xvg 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

rmsd(){
    printf "\t\tCalculating RMS Deviation........"
    if [ ! -f rmsd/$MOLEC.crystal.xvg ] ; then 
        create_dir rmsd
        cd rmsd

        if [ ! -f crystal.ndx ] ; then 
            echo '4 && ri 3-230' > selection.dat 
            echo 'q ' >> selection.dat 
            
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.tpr -o crystal.ndx >> $logFile 2>> $errFile 
            check crystal.ndx 
            fi  
        if [ ! -f $MOLEC.crystal.xvg ] ; then 
            echo 'Backbone_&_r_3-230 Backbone_&_r_3-230' | gmx rms \
                -s ../Production/$MOLEC.production.tpr \
                -f ../Production/$MOLEC.production.nopbc.xtc \
                -o $MOLEC.crystal.xvg \
                -n crystal.ndx >> $logFile 2>> $errFile 
            fi  

        check $MOLEC.crystal.xvg 
        clean 
        printf "Success\n"
        cd ../ 
    else 
        printf "Skipped\n"
        fi  
}

r_gyrate(){
    printf "\t\tCalculating radius of gyration..."
    if [ ! -f gyrate/$MOLEC.gyrate.xvg ] ; then 
        create_dir gyrate   
        cd gyrate
        if [ ! -f crystal.ndx ] ; then 
            echo '4 && ri 3-230' > selection.dat 
            echo 'q' >> selection.dat 

            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.tpr -o crystal.ndx >> $logFile 2>> $errFile 
            check crystal.ndx 
            fi 
        if [ ! -f $MOLEC.gyrate.xvg ] ; then 
            echo 'Backbone_&_r_3-230' | gmx gyrate -s ../Production/$MOLEC.production.tpr -f ../Production/$MOLEC.production.xtc -o $MOLEC.gyrate.xvg -n crystal.ndx >> $logFile 2>> $errFile 
            fi 

        check $MOLEC.gyrate.xvg 
        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
}

chi1_his148(){
    printf "\t\tCalculation chi1 of H148........." 
    if [[ ! -f chi1_his148/$MOLEC.angaver.xvg || ! -f chi1_his148/$MOLEC.angdist.xvg ]] ; then 
        create_dir chi1_his148 
        cd chi1_his148 
    
        N=`grep " N " ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CA=`grep " CA" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CB=`grep " CB" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CG=`grep " CG" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        echo "[ chi1_his148 ] " > chi1_his148.ndx 
        echo "$N   $CA   $CB   $CG   " >> chi1_his148.ndx 
        echo " " >> chi1_his148.ndx 
        check chi1_his148.ndx 

        gmx angle -f ../Production/$MOLEC.production.nopbc.xtc -type dihedral -n chi1_his148.ndx -od $MOLEC.angdist.xvg -ov $MOLEC.angaver.xvg >> $logFile 2>> $errFile 

        check $MOLEC.angaver.xvg $MOLEC.angdist.xvg
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

chi1_cnf(){
    printf "\t\tCalculation chi1 of CNF.........." 
    if [[ ! -f chi1_cnf/$MOLEC.angaver.xvg || ! -f chi1_cnf/$MOLEC.angdist.xvg ]] ; then 
        create_dir chi1_cnf 
        cd chi1_cnf 
    
        N=`grep " N " ../Production/$MOLEC.production.nopbc.gro | grep "CNF" | awk '{print$3}'`
        CA=`grep " CA" ../Production/$MOLEC.production.nopbc.gro | grep "CNF" | awk '{print$3}'`
        CB=`grep " CB" ../Production/$MOLEC.production.nopbc.gro | grep "CNF" | awk '{print$3}'`
        CG=`grep " CG" ../Production/$MOLEC.production.nopbc.gro | grep "CNF" | awk '{print$3}'`
        echo "[ chi1_cnf ] " > chi1_cnf.ndx 
        echo "$N   $CA   $CB   $CG   " >> chi1_cnf.ndx 
        echo " " >> chi1_cnf.ndx 
        check chi1_cnf.ndx 

        gmx angle -f ../Production/$MOLEC.production.nopbc.xtc -type dihedral -n chi1_cnf.ndx -od $MOLEC.angdist.xvg -ov $MOLEC.angaver.xvg >> $logFile 2>> $errFile 

        check $MOLEC.angaver.xvg $MOLEC.angdist.xvg
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

chi(){
    printf "\t\tCalculating dihedrals............" 
    if [ ! -f chi/order.xvg ] ; then 
        create_dir chi 
        cd chi 
        clean 

        gmx chi -f ../Production/$MOLEC.production.xtc \
            -s ../Production/$MOLEC.production.tpr \
            -norad \
            -rama >> $logFile 2>> $errFile 
        check order.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
protein_min
solvate
solvent_min
production_run 
if grep -sq CNF $MOLEC.pdb ; then 
    analyze_hbond
    analyze_hbond_nit
    force_calc
    force_calc_APBS
    min_dist
    sasa
    fi 
rmsd
r_gyrate
chi
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
