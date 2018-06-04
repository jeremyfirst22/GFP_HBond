#!/bin/bash

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} > [ temperature: 310K, 290K. Default = 300K]"
    exit 
}

if [ -z $1 ] ; then 
    usage 
    fi 

if [ ! -z $2 ] ; then 
    if [ $2 == 310 ] ; then 
        temp=true
        temperature=310
    elif [ $2 == 290 ] ; then 
        temp=true
        temperature=290
    else 
        echo "ERROR: $2 not a implemented temperature yet." ; exit 
    fi 
else 
    temp=false
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
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/. ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log 
errFile=$TOP/$MOLEC/$MOLEC.err 
FF=$TOP/GMXFF
forceField=amber03
FORCE_TOOLS=/Users/jfirst/force_calc_tools
if [ ! -d $FF/$forceField.ff ] ; then 
    echo ; echo "ERROR: FF not found" 
    exit
    fi  

check(){
    for arg in $@ ; do  
         if [ ! -s $arg ] ; then 
             echo ; echo "ERROR: $arg missing. Exitting" 
             exit 
             fi  
         done 
}

clean(){
    if [ -d $forceField.ff ] ; then rm -r $forceField.ff *.dat ; fi  
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi  

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi  
                                                        
    if [ ! -d $dirName/$forceField.ff ] ; then 
        if [ -d $FF/$forceField.ff ] ; then 
            cp -r $FF/$forceField.ff $dirName
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found" 
            exit 
            fi  
        fi  
}

protein_steep(){
    printf "\t\tProtein steep............................." 
    if [ ! -f Protein_steep/protein_steep.gro ] ; then 
        create_dir Protein_steep
        
        cp $MOLEC.pdb Protein_steep/.
        cd Protein_steep

        gmx pdb2gmx -f $MOLEC.pdb \
            -p $MOLEC.top \
            -ff $forceField \
            -water tip3p \
            -o $MOLEC.gro >> $logFile 2>> $errFile 
        check $MOLEC.gro 

        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -d 1.5 \
            -bt dodecahedron \
            -o boxed.gro >> $logFile 2>> $errFile
        check boxed.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c boxed.gro \
            -p $MOLEC.top \
            -o protein_steep.tpr >> $logFile 2>> $errFile 
        check protein_steep.tpr 

        gmx mdrun -deffnm protein_steep \
            -nt 128 >> $logFile 2>> $errFile 
        check protein_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi
} 

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        #cp Protein_steep/$MOLEC*.itp Solvate/. 
        cd Solvate

        gmx solvate -cp protein_steep.gro \
            -p $MOLEC.top \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        gmx grompp -f $MDP/vac_md.mdp \
            -p $MOLEC.top \
            -c solvated.gro \
            -o genion.tpr >> $logFile 2>> $errFile 
        check genion.tpr
        
        echo 'SOL' | gmx genion -s genion.tpr \
            -neutral \
            -nname 'CL' \
            -pname 'NA' \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        gmx pdb2gmx -f neutral.gro \
            -ff $forceField \
            -water tip3p \
            -p neutral.top \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        sed 's/POSRES/POSRES_IONS/' neutral_Ion2.itp > temp.itp 
        mv temp.itp neutral_Ion2.itp 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

solvent_steep(){
    printf "\t\tSolvent steep............................." 
    if [ ! -f Solvent_steep/solvent_steep.gro ] ; then 
        create_dir Solvent_steep
        
        cp Solvate/neutral.gro Solvent_steep/. 
        cp Solvate/neutral.top Solvent_steep/. 
        cp Solvate/neutral*.itp Solvent_steep/. 
        cp Solvate/posre**.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep.tpr >> $logFile 2>> $errFile 
        check solvent_steep.tpr 

        gmx mdrun -deffnm solvent_steep \
            -nt 128 >> $logFile 2>> $errFile 
        check solvent_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt(){
    printf "\t\tSolvent NVT relaxation...................." 
    if [ ! -f Solvent_nvt/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt/. 
        cp Solvent_steep/neutral.top Solvent_nvt/. 
        cp Solvent_steep/*.itp Solvent_nvt/. 
        cd Solvent_nvt

        gmx grompp -f $MDP/solvent_nvt_relax.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt \
            -nt 128 >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt_temp(){
    printf "\t\tSolvent NVT relaxation at %3i K..........." $temperature
    if [ ! -f Solvent_nvt_$temperature/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt_$temperature
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt_$temperature/. 
        cp Solvent_steep/neutral.top Solvent_nvt_$temperature/. 
        cp Solvent_steep/*.itp Solvent_nvt_$temperature/. 
        cd Solvent_nvt_$temperature

        gmx grompp -f $MDP/solvent_nvt_relax_$temperature.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt \
            -nt 128 >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt(){
    printf "\t\tSolvent NPT relaxation...................." 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_npt_relax.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt \
            -nt 128 >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt_temp(){
    printf "\t\tSolvent NPT relaxation at %3i K..........." $temperature 
    if [ ! -f Solvent_npt_$temperature/solvent_npt.gro ] ; then 
        create_dir Solvent_npt_$temperature
        
        cp Solvent_nvt_$temperature/solvent_nvt.gro Solvent_npt_$temperature/. 
        cp Solvent_nvt_$temperature/neutral.top Solvent_npt_$temperature/. 
        cp Solvent_nvt_$temperature/*.itp Solvent_npt_$temperature/. 
        cd Solvent_npt_$temperature

        gmx grompp -f $MDP/solvent_npt_relax_$temperature.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt \
            -nt 128 >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

production(){
    printf "\t\tProduction run............................" 
    if [ ! -f Production/$MOLEC.nopbc.gro ] ; then 
        create_dir Production
        
        cp Solvent_npt/neutral.top Production/.
        cp Solvent_npt/solvent_npt.gro Production/.
        cp Solvent_npt/*.itp Production/. 
        cd Production

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production.mdp \
                    -p neutral.top \
                    -c solvent_npt.gro \
                    -o $MOLEC.tpr >> $logFile 2>> $errFile 
                fi 
                check $MOLEC.tpr 

            if [ -f $MOLEC.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC \
                    -cpi $MOLEC.cpt \
                    -nt 128 >> $logFile 2>> $errFile  
            else 
                gmx mdrun -deffnm $MOLEC \
                    -nt 128 >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.gro 

        if [ ! -f $MOLEC.nopbc.xtc ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.xtc \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.gro \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

production_temp(){
    printf "\t\tProduction run at %3i K..................." $temperature
    if [ ! -f Production_$temperature/$MOLEC.nopbc.gro ] ; then 
        create_dir Production_$temperature
        
        cp Solvent_npt_$temperature/neutral.top Production_$temperature/.
        cp Solvent_npt_$temperature/solvent_npt.gro Production_$temperature/.
        cp Solvent_npt_$temperature/*.itp Production_$temperature/. 
        cd Production_$temperature

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production_$temperature.mdp \
                    -p neutral.top \
                    -c solvent_npt.gro \
                    -o $MOLEC.tpr >> $logFile 2>> $errFile 
                fi 
                check $MOLEC.tpr 

            if [ -f $MOLEC.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC \
                    -cpi $MOLEC.cpt \
                    -nt 128 >> $logFile 2>> $errFile  
            else 
                gmx mdrun -deffnm $MOLEC \
                    -nt 128 >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.gro 

        if [ ! -f $MOLEC.nopbc.xtc ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.xtc \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.gro \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond(){
    printf "\t\tAnalyzing hydrogen bonds.................." 
    if [ ! -f hbond/geometry.xvg ] ; then 
        create_dir hbond
        cp Production/solvent_npt.gro hbond/. 
        cp Production/neutral.top hbond/. 
        cp Production/*.itp hbond/. 
        cd hbond
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -op persistent.xvg \
            -or geometry.xvg \
            -oa hb_count.xvg \
            -onwr nw_geometry.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg geometry.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 
   
#        echo "r CNF & a NH" > selection.dat 
#        echo "r SOL" >> selection.dat 
#        echo "!r CNF" >> selection.dat  
#        echo "q" >> selection.dat  
#
#        touch empty.ndx 
#        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
#            -n empty.ndx \
#            -o index.ndx >> $logFile 2>> $errFile 
#        check index.ndx 
#
#        echo '0 1 0' | gmx hbond -f ../Production/$MOLEC.xtc \
#            -s ../Production/$MOLEC.tpr \
#            -n index.ndx \
#            -shell 1.0 \
#            -dist wat_hbdist.xvg \
#            -ang wat_hbang.xvg \
#            -num wat_hbnum.xvg >> $logFile 2>> $errFile 
#            #-ac wat_hbac.xvg \
#            #-life wat_hblife.xvg \
#        check wat_hbnum.xvg 
#
#        echo '0 2 0' | gmx hbond -f ../Production/$MOLEC.xtc \
#            -s ../Production/$MOLEC.tpr \
#            -n index.ndx \
#            -shell 1.0 \
#            -ac all_hbac.xvg \
#            -dist all_hbdist.xvg \
#            -ang all_hbang.xvg \
#            -life all_hblife.xvg \
#            -num all_hbnum.xvg >> $logFile 2>> $errFile 
#        check all_hbnum.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_with_ca(){
    printf "\t\tAnalyzing hydrogen bonds.................." 
    if [ ! -f hbond_with_ca/geometry.xvg ] ; then 
        create_dir hbond_with_ca
        cp Production/solvent_npt.gro hbond_with_ca/. 
        cp Production/neutral.top hbond_with_ca/. 
        cp Production/*.itp hbond_with_ca/. 
        cd hbond_with_ca
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -op persistent.xvg \
            -or geometry.xvg \
            -oa hb_count.xvg \
            -onwr nw_geometry.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg geometry.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_with_ca_temp(){
    printf "\t\tAnalyzing hydrogen bonds.................." 
    if [ ! -f hbond_with_ca_$temperature/geometry.xvg ] ; then 
        create_dir hbond_with_ca_$temperature
        cp Production_$temperature/solvent_npt.gro hbond_with_ca_$temperature/. 
        cp Production_$temperature/neutral.top hbond_with_ca_$temperature/. 
        cp Production_$temperature/*.itp hbond_with_ca_$temperature/. 
        cd hbond_with_ca_$temperature
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production_$temperature/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -op persistent.xvg \
            -or geometry.xvg \
            -oa hb_count.xvg \
            -onwr nw_geometry.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg geometry.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

fit_hbond_with_ca(){
    printf "\t\tFitting hbond analysia...................." 
    if [ ! -f hbond_with_ca/geometry.xvg ] ; then 
        printf "Skipping\n"  
        break 
        fi 

    if [ ! -f ~/normal_distribution/tiltAngle ] ; then 
        printf "Skipping\n" 
        break 
        fi     

    if [ ! -f fit_hbond_with_ca/theta1.poly ] ; then 
        create_dir fit_hbond_with_ca
        cp hbond_with_ca/*geometry.xvg fit_hbond_with_ca/. 
        cd fit_hbond_with_ca
        clean 

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > clean_dist.dat 
        check clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_dist.dat -o dist.his -g dist.gaus -p dist.poly 
        check dist.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > clean_theta1.dat 
        check clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta1.dat -o theta1.his -g theta1.gaus -p theta1.poly -n 81 
        check theta1.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > clean_theta2.dat 
        check clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta2.dat -o theta2.his -g theta2.gaus -p theta2.poly 
        check theta2.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > nw_clean_dist.dat 
        check nw_clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_dist.dat -o nw_dist.his -g nw_dist.gaus -p nw_dist.poly 
        check nw_dist.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > nw_clean_theta1.dat 
        check nw_clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta1.dat -o nw_theta1.his -g nw_theta1.gaus -p nw_theta1.poly 
        check nw_theta1.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > nw_clean_theta2.dat 
        check nw_clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta2.dat -o nw_theta2.his -g nw_theta2.gaus -p nw_theta2.poly 
        check nw_theta2.poly

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

fit_hbond_with_ca_temp(){
    printf "\t\tFitting hbond analysia...................." 
    if [ ! -f hbond_with_ca_$temperature/geometry.xvg ] ; then 
        printf "Skipping\n"  
        break 
        fi 

    if [ ! -f ~/normal_distribution/tiltAngle ] ; then 
        printf "Skipping\n" 
        break 
        fi     

    if [ ! -f fit_hbond_with_ca_$temperature/theta1.poly ] ; then 
        create_dir fit_hbond_with_ca_$temperature
        cp hbond_with_ca_$temperature/*geometry.xvg fit_hbond_with_ca_$temperature/. 
        cd fit_hbond_with_ca_$temperature
        clean 

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > clean_dist.dat 
        check clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_dist.dat -o dist.his -g dist.gaus -p dist.poly 
        check dist.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > clean_theta1.dat 
        check clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta1.dat -o theta1.his -g theta1.gaus -p theta1.poly -n 81 
        check theta1.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > clean_theta2.dat 
        check clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta2.dat -o theta2.his -g theta2.gaus -p theta2.poly 
        check theta2.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > nw_clean_dist.dat 
        check nw_clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_dist.dat -o nw_dist.his -g nw_dist.gaus -p nw_dist.poly 
        check nw_dist.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > nw_clean_theta1.dat 
        check nw_clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta1.dat -o nw_theta1.his -g nw_theta1.gaus -p nw_theta1.poly 
        check nw_theta1.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > nw_clean_theta2.dat 
        check nw_clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta2.dat -o nw_theta2.his -g nw_theta2.gaus -p nw_theta2.poly 
        check nw_theta2.poly

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_1(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 1..." 
    if [ ! -f hbond_1/persistent.xvg ] ; then 
        create_dir hbond_1
        cp Production/solvent_npt.gro hbond_1/. 
        cp Production/neutral.top hbond_1/. 
        cp Production/*.itp hbond_1/. 
        cd hbond_1
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 1 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_2(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 2..." 
    if [ ! -f hbond_2/persistent.xvg ] ; then 
        create_dir hbond_2
        cp Production/solvent_npt.gro hbond_2/. 
        cp Production/neutral.top hbond_2/. 
        cp Production/*.itp hbond_2/. 
        cd hbond_2
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 2 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_2_temp(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 2..." 
    if [ ! -f hbond_2_$temperature/persistent.xvg ] ; then 
        create_dir hbond_2_$temperature
        cp Production_$temperature/solvent_npt.gro hbond_2_$temperature/. 
        cp Production_$temperature/neutral.top hbond_2_$temperature/. 
        cp Production_$temperature/*.itp hbond_2_$temperature/. 
        cd hbond_2_$temperature
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production_$temperature/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 2 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_3(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 3..." 
    if [ ! -f hbond_3/persistent.xvg ] ; then 
        create_dir hbond_3
        cp Production/solvent_npt.gro hbond_3/. 
        cp Production/neutral.top hbond_3/. 
        cp Production/*.itp hbond_3/. 
        cd hbond_3
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 3 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_4(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 4..." 
    if [ ! -f hbond_4/persistent.xvg ] ; then 
        create_dir hbond_4
        cp Production/solvent_npt.gro hbond_4/. 
        cp Production/neutral.top hbond_4/. 
        cp Production/*.itp hbond_4/. 
        cd hbond_4
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 4 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_5(){
    printf "\t\tAnalyzing hydrogen bonds with forgive 5..." 
    if [ ! -f hbond_5/persistent.xvg ] ; then 
        create_dir hbond_5
        cp Production/solvent_npt.gro hbond_5/. 
        cp Production/neutral.top hbond_5/. 
        cp Production/*.itp hbond_5/. 
        cd hbond_5
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNF solvent_npt.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF solvent_npt.gro | grep NH | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -b 10000 \
            -select 'not resname CNF and (same residue as within 0.5 of resname CNF and name NH)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 5 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

fit_hbond(){
    printf "\t\tFitting hbond analysia...................." 
    if [ ! -f hbond/geometry.xvg ] ; then 
        printf "Skipping\n"  
        break 
        fi 

    if [ ! -f ~/normal_distribution/tiltAngle ] ; then 
        printf "Skipping\n" 
        break 
        fi     

    if [ ! -f fit_hbond/nw_theta2.poly ] ; then 
        create_dir fit_hbond
        cp hbond/*geometry.xvg fit_hbond/. 
        cd fit_hbond
        clean 

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > clean_dist.dat 
        check clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_dist.dat -o dist.his -g dist.gaus -p dist.poly 
        check dist.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > clean_theta1.dat 
        check clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta1.dat -o theta1.his -g theta1.gaus -p theta1.poly -n 81 
        check theta1.poly

        check geometry.xvg 
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > clean_theta2.dat 
        check clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta2.dat -o theta2.his -g theta2.gaus -p theta2.poly 
        check theta2.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > nw_clean_dist.dat 
        check nw_clean_dist.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_dist.dat -o nw_dist.his -g nw_dist.gaus -p nw_dist.poly 
        check nw_dist.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > nw_clean_theta1.dat 
        check nw_clean_theta1.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta1.dat -o nw_theta1.his -g nw_theta1.gaus -p nw_theta1.poly 
        check nw_theta1.poly

        check nw_geometry.xvg 
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > nw_clean_theta2.dat 
        check nw_clean_theta2.dat 
       
        /Users/jfirst/normal_distribution/tiltAngle -f nw_clean_theta2.dat -o nw_theta2.his -g nw_theta2.gaus -p nw_theta2.poly 
        check nw_theta2.poly

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

sasa(){
    printf "\t\tAnalyzing SASA of CNF....................." 
    if [ ! -f sasa/aromatic.xvg ] ; then 
        create_dir sasa
        cd sasa

        if [ ! -s cnf_area.xvg ] ; then 
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNF' \
                -ndots 240 \
                -or cnf_resarea.xvg \
                -o cnf_area.xvg >> $logFile 2>> $errFile 
            fi 
        check cnf_area.xvg 

        if [ ! -s nh_ct_area.xvg ] ; then 
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNF and (name NH or name CT)' \
                -ndots 240 \
                -o nh_ct_area.xvg >> $logFile 2>> $errFile
            fi 
        check nh_ct_area.xvg 

        if [ ! -s nit_4_atoms_area.xvg ] ; then 
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNF and (name NH or name CT or name CZ or name CE1 or name CE2)' \
                -ndots 240 \
                -o nit_4_atoms_area.xvg >> $logFile 2>> $errFile
            fi 
        check nit_4_atoms_area.xvg 

        if [ ! -s sidechain.xvg ] ; then 
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'group "SideChain" and resname CNF' \
                -ndots 240 \
                -o sidechain.xvg >> $logFile 2>> $errFile
            fi 
        check sidechain.xvg 

        if [ ! -s aromatic.xvg ] ; then 
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNF and not name N H C O CA HA CB CB1 CB2' \
                -ndots 240 \
                -o aromatic.xvg >> $logFile 2>> $errFile
            fi 
        check aromatic.xvg 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi 
}

sorient(){
    printf "\t\tCalculating orientation of solvent........" 
    if [ ! -f sorient/sori.xvg ] ; then 
        create_dir sorient
        cd sorient
        clean 
        
        touch empty.ndx 
        echo "r CNF & a NH" > selection.dat 
        echo "r SOL" >> selection.dat 
        echo "q" >> selection.dat 
        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 1' | gmx sorient -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -no snor.xvg \
            -ro sord.xvg \
            -co scum.xvg \
            -rc scount.xvg \
            -o sori.xvg >> $logFile 2>> $errFile 
        check sori.xvg snor.xvg sord.xvg scum.xvg scount.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi 
}

rdf(){
    printf "\t\tCalculating radial distribution funcs....." 
    if [ ! -f rdf/wat_cnf.xvg ] ; then 
        create_dir rdf
        cd rdf
        clean 

        touch empty.ndx 
        echo "r CNF & a NH" > selection.dat 
        echo "r SOL & a OW" >> selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o wat_cnf.xvg >> $logFile 2>> $errFile 
        check wat_cnf.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi 

}

force_calc(){
    printf "\n\t\tCalculating force:\n" 
    if [[ ! -f force_calc/$MOLEC.solvent_rxn_field.projected.xvg || ! -f force_calc/$MOLEC.external_field.projected.xvg || ! -f force_calc/$MOLEC.total_field.xvg ]] ; then 

        if [ ! -f $FORCE_TOOLS/g_insert_dummy_atom ] ; then 
            printf "\t\t\tERROR: Force tools not found. Skipping force calc\n" 
            return  
            fi 

        create_dir force_calc
        cp Production/*.itp force_calc/. 
        cp Production/neutral.top force_calc/. 
        cp Production/solvent_npt.gro force_calc/. 
        cp Solvate/neutral.gro force_calc/. 
        cd force_calc 

        ## We use version 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
        ## We allow for warnings, since we are generated .tpr from a gromacs 5 mdp file. We are only inserting
        ## atoms this should not matter. 
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 
        
        CT=`grep CNF ../Production/$MOLEC.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF ../Production/$MOLEC.gro | grep NH | awk '{print $3}'`
        #echo $CT $NH
        
        printf "\t\t\tInserting dummy atoms............................" 
        if [ ! -s $MOLEC.with_dummy.xtc ] ; then 
            source /usr/local/gromacs/bin/GMXRC
            $FORCE_TOOLS/g_insert_dummy_atom -f ../Production/$MOLEC.nopbc.xtc \
                -s v4.tpr \
                -a1 $CT \
                -a2 $NH \
                -o $MOLEC.with_dummy.xtc >> $logFile 2>> $errFile 
            check $MOLEC.with_dummy.xtc

        ## We use the initial configuration so that titration states are conserved (ie, at the end of the production run, pdb2gmx might assign a different titration state to a histidine, which causes it to fail. 
            if [ ! -s $MOLEC.with_dummy.gro ] ; then 
                echo 'Protein System' | gmx trjconv -f neutral.gro \
                    -s v4.tpr \
                    -center \
                    -ur compact \
                    -pbc mol \
                    -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
                check $MOLEC.nopbc.gro 

                source /usr/local/gromacs/bin/GMXRC
                $FORCE_TOOLS/g_insert_dummy_atom -f $MOLEC.nopbc.gro \
                    -s v4.tpr \
                    -a1 $CT \
                    -a2 $NH \
                    -o $MOLEC.with_dummy.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.with_dummy.gro 

            if [ ! -s $MOLEC.with_dummy.top ] ; then 
                gmx pdb2gmx -f $MOLEC.with_dummy.gro \
                    -water tip3p \
                    -ff amber03 \
                    -p $MOLEC.with_dummy.top \
                    -o $MOLEC.with_dummy.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.with_dummy.top 
            printf "Done\n" 
        else 
            printf "Skipped\n" 
            fi 
        
        ##Find new atom numbers 
        CT=`grep CNF $MOLEC.with_dummy.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF $MOLEC.with_dummy.gro | grep NH | awk '{print $3}'`

        echo "[ probe ]" > probe.ndx 
        echo "$CT $NH" >> probe.ndx 

        echo "[ protein ]" > protein.ndx 
        grep -v TCHG $MOLEC.with_dummy.gro | grep -v SOL | grep -v HOH | grep -v NA | grep -v CL | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 

        cp $MOLEC.with_dummy.top $MOLEC.total_field.top 

        if [ ! -s $MOLEC.solvent_rxn_field.top ] ; then 
            $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top protein.ndx $MOLEC.solvent_rxn_field.top >> $logFile 2>> $errFile 
            fi 

        if [ ! -s $MOLEC.external_field.top ] ; then 
            $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top probe.ndx $MOLEC.external_field.top >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.total_field.top $MOLEC.external_field.top $MOLEC.solvent_rxn_field.top 

        for field in total_field external_field solvent_rxn_field ; do 
            printf "\t\t%20s..." $field 

            ##Extract forces 
            if [ ! -s $MOLEC.$field.projected.xvg ] ; then 
                printf "forces..." 
                if [ ! -s $MOLEC.$field.xvg ] ; then 
                    if [ ! -s $MOLEC.$field.tpr ] ; then 
                        gmx grompp -f $MDP/rerun.mdp \
                            -p $MOLEC.$field.top \
                            -c $MOLEC.with_dummy.gro \
                            -o $MOLEC.$field.tpr  >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$field.tpr 
 
                    if [ ! -s $MOLEC.$field.trr ] ; then 
                        gmx mdrun -s $MOLEC.$field.tpr \
                            -rerun $MOLEC.with_dummy.xtc \
                            -deffnm $MOLEC.$field >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$field.trr 

                    echo 2 | gmx traj -f $MOLEC.$field.trr \
                        -s $MOLEC.$field.tpr \
                        -xvg none \
                        -of $MOLEC.$field.xvg >> $logFile 2>> $errFile 
                    rm $MOLEC.$field.trr 
                fi 
                check $MOLEC.$field.xvg 

                ##extract postions for bond vector
                printf "positions..." 
                if [ ! -s $MOLEC.positions.xvg ] ; then 
                    gmx traj -f $MOLEC.with_dummy.xtc \
                        -s $MOLEC.$field.tpr \
                        -n probe.ndx \
                        -xvg none \
                        -ox $MOLEC.positions.xvg >> $logFile 2>> $errFile 
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
        cd ../
    else 
        printf "\t\t\t\t  ............Skipped\n" 
        fi 
    printf "\n" 
}

minimage(){
    printf "\t\tCalculating minimum image................." 
    if [ ! -f minimage/mindist.xvg ] ; then 
        create_dir minimage
        cd minimage
        clean 
        
        echo 'Protein' | gmx mindist -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pi \
            -od mindist.xvg >> $logFile 2>> $errFile 
            check mindist.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

mindist(){
    printf "\t\tCalculating mindist to water.............." 
    if [ ! -f mindist/mindist.xvg ] ; then 
        create_dir mindist
        cd mindist
        clean 
        
        touch empty.ndx 
        echo "r CNF & a NH" > selection.dat 
        echo "r SOL & a OW" >> selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.gro \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 1' | gmx mindist -f ../Production/$MOLEC.nopbc.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -od mindist.xvg >> $logFile 2>> $errFile 
        check mindist.xvg

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

rmsd(){
    printf "\t\tCalculating RMSD.........................." 
    if [ ! -f rmsd/without_ter.xvg ] ; then 
        create_dir rmsd
        cd rmsd
        clean 

        echo 'Backbone Backbone' | gmx rms -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -b 10000 \
            -o backbone.xvg >> $logFile 2>> $errFile 
        check backbone.xvg 

        echo '"Backbone" & ri 11-225' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_11_225 Backbone_&_r_11_225" | gmx rms -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -b 10000 \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
}

rmsd_temp(){
    printf "\t\tCalculating RMSD.........................." 
    if [ ! -f rmsd_$temperature/without_ter.xvg ] ; then 
        create_dir rmsd_$temperature
        cd rmsd_$temperature
        clean 

        echo 'Backbone Backbone' | gmx rms -f ../Production_$temperature/$MOLEC.xtc \
            -s ../Production_$temperature/$MOLEC.tpr \
            -b 10000 \
            -o backbone.xvg >> $logFile 2>> $errFile 
        check backbone.xvg 

        echo '"Backbone" & ri 11-225' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production_$temperature/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_11_225 Backbone_&_r_11_225" | gmx rms -f ../Production_$temperature/$MOLEC.xtc \
            -s ../Production_$temperature/$MOLEC.tpr \
            -n index.ndx \
            -b 10000 \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
}

chi(){
    printf "\t\tCalculating phi psi chi dihedrals........." 
    if [ ! -f chi/order.xvg ] ; then 
        create_dir chi
        cd chi
        clean 

        gmx chi -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -b 10000 \
            -norad \
            -rama >> $logFile 2>> $errFile 
        check order.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

chi_temp(){
    printf "\t\tCalculating phi psi chi dihedrals........." 
    if [ ! -f chi_$temperature/order.xvg ] ; then 
        create_dir chi_$temperature
        cd chi_$temperature
        clean 

        gmx chi -f ../Production_$temperature/$MOLEC.xtc \
            -s ../Production_$temperature/$MOLEC.tpr \
            -b 10000 \
            -norad \
            -rama >> $logFile 2>> $errFile 
        check order.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

rgyr(){
    printf "\t\tCalculating radius of gyration............" 
    if [ ! -f rgyr/without_ter.xvg ] ; then 
        create_dir rgyr
        cd rgyr
        clean 

        echo 'Protein' | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o gyrate.xvg >> $logFile 2>> $errFile 
        check gyrate.xvg

        echo '"Backbone" & ri 11-225' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_11_225 Backbone_&_r_11_225" | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -b 10000 \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

rgyr_temp(){
    printf "\t\tCalculating radius of gyration............" 
    if [ ! -f rgyr_$temperature/without_ter.xvg ] ; then 
        create_dir rgyr_$temperature
        cd rgyr_$temperature
        clean 

        echo 'Protein' | gmx gyrate -f ../Production_$temperature/$MOLEC.xtc \
            -s ../Production_$temperature/$MOLEC.tpr \
            -o gyrate.xvg >> $logFile 2>> $errFile 
        check gyrate.xvg

        echo '"Backbone" & ri 11-225' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production_$temperature/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_11_225 Backbone_&_r_11_225" | gmx gyrate -f ../Production_$temperature/$MOLEC.xtc \
            -s ../Production_$temperature/$MOLEC.tpr \
            -n index.ndx \
            -b 10000 \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
protein_steep
solvate
solvent_steep
if ! $temp ; then 
    solvent_nvt
    solvent_npt
    production 
    if grep -sq CNF $MOLEC.pdb ; then 
    #    hbond 
        hbond_with_ca
    #    hbond_1
        hbond_2
    #    hbond_3
    #    hbond_4
    #    hbond_5
        fit_hbond_with_ca
       sasa
    #    mindist 
    #    sorient
    #    rdf
    #    force_calc
        chi
        fi 
    rmsd 
    rgyr
else 
    solvent_nvt_temp
    solvent_npt_temp
    production_temp
    if grep -sq CNF $MOLEC.pdb ; then 
        hbond_with_ca_temp
        hbond_2_temp
        fit_hbond_with_ca_temp
        fi 
    chi_temp
    rmsd_temp
    rgyr_temp
fi 
#minimage
#rgyr
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
