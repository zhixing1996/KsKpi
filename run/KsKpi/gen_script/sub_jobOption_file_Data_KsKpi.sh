#/bin/bash

BOSS=$1

if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

for ECM in ${ECMS[@]}; do
    cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/jobs_text/data/$ECM
    find . -name "*.bosslog" | xargs rm
    find . -name "*.bosserr" | xargs rm
    NUM_UP=$(ls -l | grep "txt" | wc -l)
    boss.condor -g physics -n $NUM_UP jobOptions_KsKpi_data_$ECM-%{ProcId}.txt
done
