#/bin/bash

BOSS=$1

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

for ECM in ${ECMS[@]}; do
    cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/rootfile/data/$ECM
    if [ "$ECM" == "3097" ]; then
        rm -rf KsKpi_data_3097-0.root
    fi
    rm -rf KsKpi_data_${ECM}.root
    hadd KsKpi_data_${ECM}.root *.root
done
