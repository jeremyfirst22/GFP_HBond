#!/bin/bash


for file in `ls *.pdb` ; do 
    echo $file 

    if grep -sq "GLH   221" $file ; then 
        echo "GLH already found! Skipping" 
    else  
        print "Fixing!"
        sed "s/GLU   221/GLH   221/" $file > new_$file
        sed "s/HIS    76/HIE    76/" new_$file > new_2_$file
        sed "s/HIS    80/HIE    80/" new_2_$file > new_3_$file
        sed "s/HIS   138/HIE   138/" new_3_$file > new_4_$file
        sed "s/HIS   147/HID   147/" new_4_$file > new_5_$file
        sed "s/HIS   168/HIE   168/" new_5_$file > new_6_$file
        sed "s/HIS   180/HID   180/" new_6_$file > new_7_$file
        sed "s/HIS   198/HID   198/" new_7_$file > new_8_$file
        sed "s/HIS   216/HID   216/" new_8_$file > new_9_$file
        sed "s/HIS   230/HIS   230/" new_9_$file > new_10_$file
        
        mv new_10_$file $file 
        rm -v new_*.pdb 
        fi 

    grep "GLH " $file | awk '{print $4"\t"$5}' | uniq 
    grep "HI[SDIE] " $file | awk '{print $4"\t"$5}' | uniq 


    done 

