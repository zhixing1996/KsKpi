#/bin/bash

BOSS=$1

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
    PATH="/besfs4/offline/data/705-1/jpsi/round02"
fi

FILENAME="Gen_FileList_Data_"$BOSS
echo "#!/usr/bin/env bash" > $FILENAME
echo "cd /besfs/groups/cal/dedx/$USER/bes/KsKpi" >> $FILENAME 
for ECM in ${ECMS[@]}; do
    echo "rm -r /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/samples/data/$ECM/*txt" >> $FILENAME
    echo "./python/get_samples.py $PATH/dst ./run/KsKpi/samples/data/$ECM/data_${ECM}_KsKpi.txt 20G" >> $FILENAME
done
