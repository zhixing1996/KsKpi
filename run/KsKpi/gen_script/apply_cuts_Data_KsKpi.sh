#/bin/bash

BOSS=$1
SOURCE_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/rootfile/data
ANA_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/anaroot/data
LOG_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/logfile
mkdir -p $ANA_PATH
mkdir -p $LOG_PATH

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

FILENAME="Apply_Cuts_Data_"$BOSS
echo "#!/usr/bin/env bash" > $FILENAME
echo "cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/python" >> $FILENAME 
for ECM in ${ECMS[@]}; do
    mkdir -p $ANA_PATH/$ECM
    echo "./apply_cuts.py $SOURCE_PATH/$ECM/KsKpi_data_${ECM}.root $ANA_PATH/$ECM/data_${ECM}_signal.root" >> $FILENAME
done
