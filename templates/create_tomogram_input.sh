#!/bin/bash

## DEFINITIONS ##

find_index()
{
    echo $1  | grep -bo $2 | sed 's/:.*$//'
}


#####  Input params  #####

# Parent folder containing ETOMO projects folder
projectData='tomograms'
prefix='TS_'
tomoDigits=2

# Output tomogram input data starfile
tomoParamStar="tomograms_input.star"

# Create Orderlist files from .mdoc data
add_orderlist=1

dose=2.7

# Create culled tomo files
culled_files=0
culled_folder=croppedtomos

##########################

prefixLen=${#prefix}
tomoList=""
for name in ${projectData}/TS_*; 
do
    idx=$(($(find_index $name TS_) + $prefixLen ))
    echo my idx = $idx
    num=${name:$idx}; tomoList="$tomoList $num"; 
    echo num = $num
    
done

#tomoList=$(echo $tomoList | xargs -n1 | sort -g | xargs)
echo tomoList=$tomoList



# Coordinates starfile to extract the rlnTomoName labels
tomoList=""
for name in ${projectData}/TS_*; do num=${name:13}; tomoList="$tomoList $num"; done
tomoList=$(echo $tomoList | xargs -n1 | sort -g | xargs)
echo tomoList=$tomoList
# Number of digits in rlnTomoName label



if [[ ${culled_files} -eq 1 ]]
then
    mkdir -p ${culled_folder}
fi
mkdir -p orderLists
rm -f ${tomoParamStar}
touch ${tomoParamStar}

echo > $tomoParamStar
echo "# version 3001" >> $tomoParamStar
echo >> $tomoParamStar
echo "data_global" >> $tomoParamStar
echo >> $tomoParamStar
echo loop_ >> $tomoParamStar
echo _rlnTomoName >> $tomoParamStar
echo _rlnTomoTiltSeriesName >> $tomoParamStar
echo _rlnTomoImportCtfFindFile >> $tomoParamStar
#echo _rlnTomoImportCtfPlotterFile >> $tomoParamStar
echo _rlnTomoImportImodDir >> $tomoParamStar
echo _rlnTomoImportFractionalDose >> $tomoParamStar
if [[ ${culled_files} -eq 1 ]]
then
    echo _rlnTomoImportCulledFile >> $tomoParamStar
fi  
if [[ ${add_orderlist} -eq 1 ]]
then
    echo _rlnTomoImportOrderList >> $tomoParamStar
fi

echo >> $tomoParamStar

for tomonum in $tomoList
do
#    tomonum=${tomoname: -${tomoDigits}:${tomoDigits}}
    tomoname=TS_$tomonum
    tomoFolder=$projectData/TS_$tomonum
    #dose=$(grep -R "${tomonum}" ../tomo_dose.lst  | awk -F',' '{print $2}')
    
    printf " ${tomoname}   " >> $tomoParamStar
    printf " ${tomoFolder}/TS_${tomonum}_aligned.st:mrc " >> $tomoParamStar
#    printf " ${tomoFolder}/ctfplotter/TS_${tomonum}_output.txt   " >> $tomoParamStar
    printf " ${tomoFolder}/ctffind4/TS_${tomonum}_output.txt   " >> $tomoParamStar
        printf " ${tomoFolder}   " >> $tomoParamStar
    printf " ${dose}   " >> $tomoParamStar
     
    if [[ ${culled_files} -eq 1 ]]
    then
        printf " ${culled_folder}/TS_${tomonum}.mrc " >> $tomoParamStar
    fi 
    if [[ ${add_orderlist} -eq 1 ]]
    then
        # First we create the order_list
        grep "TiltAngle" ${tomoFolder}/TS_${tomonum}.st.mdoc | awk -F= '{print NR "," sprintf("%.0f",$2)}' > orderLists/order_list_${tomonum}.csv
        printf " orderLists/order_list_${tomonum}.csv " >> $tomoParamStar
    fi
    printf "\n" >> $tomoParamStar

done

echo >> $tomoParamStar 

exit 0 
