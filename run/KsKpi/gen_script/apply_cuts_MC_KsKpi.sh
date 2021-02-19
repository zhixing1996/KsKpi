#/bin/bash

BOSS=$1
TYPE=$2
SOURCE_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/rootfile/mc/$TYPE
ANA_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/anaroot/mc/$TYPE
LOG_PATH=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/logfile
mkdir -p $ANA_PATH
mkdir -p $LOG_PATH

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

FILENAME="Apply_Cuts_${TYPE}_"$BOSS
echo "#!/usr/bin/env bash" > $FILENAME
echo "cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/python" >> $FILENAME 
for ECM in ${ECMS[@]}; do
    mkdir -p $ANA_PATH/$ECM
    echo "./apply_cuts.py $SOURCE_PATH/$ECM/KsKpi_mc_${TYPE}_${ECM}.root $ANA_PATH/$ECM/mc_${TYPE}_${ECM}_signal.root" >> $FILENAME
done
