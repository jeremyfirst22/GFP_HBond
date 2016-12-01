#!/bin/bash

molecList="GFP_Y92X GFP_F114X GFP_Y143X GFP_F145X GFP_Y149X GFP_F165X GFP_N212X GFP_M218X" 

for molec in $molecList ; do 
    if [ ! -f StartingStructures/$molec.pdb ] ; then 
        echo "$molec.pdb not found!" 
        continue 
        fi 
    
    printf "\n\n\t\t*** $molec ***\n\n" 


    echo "#!/bin/bash" > submit_$molec
    echo >> submit_$molec
    echo "#SBATCH -J $molec " >> submit_$molec
    echo "#SBATCH -o $molec.o%j" >> submit_$molec
    echo "#SBATCH -n 16 " >> submit_$molec
    echo "#SBATCH -p normal " >> submit_$molec
    echo "#SBATCH -t 48:00:00" >> submit_$molec
    echo "#SBATCH -A Ras" >> submit_$molec
    echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_$molec
    echo "#SBATCH --mail-type=all" >> submit_$molec
    
    echo >> submit_$molec
    echo "module load boost " >> submit_$molec
    echo "module load cxx11 " >> submit_$molec
    echo "module load gromacs " >> submit_$molec 
   
    echo >> submit_$molec
    echo "bash run_GFP_hbond.sh StartingStructures/${molec}.pdb" >> submit_$molec
   
    sbatch submit_$molec
#    bash run_GFP_hbond.sh StartingStructures/$molec.pdb 


    done 
