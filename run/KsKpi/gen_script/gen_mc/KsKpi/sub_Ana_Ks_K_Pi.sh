#!/bin/sh
BOSS=$1
if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
fi

WORKAREA=/besfs/groups/cal/dedx/$USER/bes/KsKpi
for ECM in ${ECMS[@]}; do
    mkdir -p $WORKAREA/run/KsKpi/rootfile/mc/KsKpi/$ECM
    mkdir -p $WORKAREA/run/KsKpi/jobs_text/mc/KsKpi/$ECM
    cd $WORKAREA/run/KsKpi/jobs_text/mc/KsKpi/$ECM
    rm -rf mc_Ks_K_Pi_HELAMP_$ECM*txt
    cp -rf $WORKAREA/python/make_mc.py ./
    cp -rf $WORKAREA/python/tools.py ./
    ./make_mc.py /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/dst/mc/KsKpi/$ECM mc Ks_K_Pi HELAMP KsKpi $ECM 2
    cp -rf /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/gen_script/gen_mc/subAna.sh ./
    rm -rf *boss*
    rm -rf $WORKAREA/run/KsKpi/rootfile/mc/KsKpi/$ECM/*root
    ./subAna.sh $ECM Ks_K_Pi_HELAMP
done
