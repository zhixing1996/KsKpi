#/bin/bash
BOSS=$1
if [ "$BOSS" = "705" ]; then
    ECMS=("3097")
    CMS=("3.0971873")
    RUNNO_LOW=("9947")
    RUNNO_UP=("10878")
fi

COUN=0
FILENAME="Sub_jobOption_"$BOSS
for ECM in ${ECMS[@]}; do
    JobText_SaveDir=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/jobs_text/mc/KsKpi/$ECM
    Script_Dir=/besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/gen_script/gen_mc/KsKpi
    mkdir -p $JobText_SaveDir
    rm -rf $JobText_SaveDir/jobOptions*txt
    rm -rf $JobText_SaveDir/subSimRec_*.sh
    rm -rf $JobText_SaveDir/fort.*
    rm -rf $JobText_SaveDir/phokhara*
    cd $JobText_SaveDir
    cp -r $Script_Dir/jobOptions_sim_sig_Ks_K_Pi_HELAMP_tempE_705.sh jobOptions_sim_sig_Ks_K_Pi_HELAMP_${ECM}_705.sh
    cp -r $Script_Dir/jobOptions_rec_sig_Ks_K_Pi_HELAMP_tempE_705.sh jobOptions_rec_sig_Ks_K_Pi_HELAMP_${ECM}_705.sh
    sh jobOptions_sim_sig_Ks_K_Pi_HELAMP_${ECM}_705.sh 0 39 $ECM ${CMS[$COUNT]} 5000 ${RUNNO_LOW[$COUNT]} ${RUNNO_UP[$COUNT]}
    sh jobOptions_rec_sig_Ks_K_Pi_HELAMP_${ECM}_705.sh 0 39 $ECM
    rm -rf /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/rtraw/KsKpi/$ECM/*.rtraw
    rm -rf /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/dst/KsKpi/$ECM/*.dst
    cp -rf /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/gen_script/gen_mc/subSimRec.sh ./
    sh subSimRec.sh jobOptions_sim_sig_Ks_K_Pi_HELAMP_$ECM jobOptions_rec_sig_Ks_K_Pi_HELAMP_$ECM subSimRec_Ks_K_Pi_$ECM 0 39
    COUNT=$((${COUNT}+1))
done
