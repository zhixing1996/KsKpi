#/bin/bash

BOSS=$1
TYPE=$2

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

for ECM in ${ECMS[@]}; do
    cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/rootfile/mc/$TYPE/$ECM
    rm -rf KsKpi_mc_${TYPE}_${ECM}.root
    hadd KsKpi_mc_${TYPE}_${ECM}.root *.root
done
