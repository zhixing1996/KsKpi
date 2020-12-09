#!/usr/bin/env bash

# Main driver to submit jobs 
# Author SHI Xin <shixin@ihep.ac.cn>
# Modified by JING Maoqiang <jingmq@ihep.ac.cn>
# Created [2019-12-11 Dec 14:56]

usage() {
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-9s  %-40s"  "0.1"       "[run on data sample for KsKpi]"
    printf "\n\t%-9s  %-40s"  "0.2"       "[generate and run on signal MC sample]"
    printf "\n\n" 
}

usage_0_1() {
    printf "\n\t%-9s  %-40s"  ""          ""   
    printf "\n\t%-9s  %-40s"  "0.1.1"     "Split data sample with each group 5G"
    printf "\n\t%-9s  %-40s"  "0.1.2"     "Generate Condor jobs on data ---- 1"
    printf "\n\t%-9s  %-40s"  "0.1.3"     "Test for data"
    printf "\n\t%-9s  %-40s"  "0.1.4"     "Submit Condor jobs on data ---- 2"
    printf "\n\t%-9s  %-40s"  "0.1.5"     "Synthesize data root files"
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_2() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.2.1"     "Generate MC samples ---- Simulation && Reconstruction"
    printf "\n\t%-9s  %-40s"  "0.2.2"     "Generate MC samples ---- Event Selection"
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

if [[ $# -eq 0 ]]; then
    usage
    echo "Please enter your option: "
    read option
else
    option=$1
fi

sub_0_1() {

case $option in
    
    # --------------------------------------------------------------------------
    #  Data  
    # --------------------------------------------------------------------------

    0.1.1) echo "Split data sample with each group 20G ..."
           cd ./run/KsKpi/gen_script
           ./make_file_list_Data_KsKpi.sh 705
           chmod u+x Gen_FileList_Data_705
           bash Gen_FileList_Data_705
           rm -r Gen_FileList_Data_705
           cd /besfs/groups/cal/dedx/$USER/bes/KsKpi
	       ;;

    0.1.2) echo "Generate Condor jobs on data ---- 1..." 
	       cd ./run/KsKpi/gen_script
	       ./make_jobOption_file_Data_KsKpi.sh 705
	       cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/jobs_text/data/3097
	       cp -r jobOptions_KsKpi_data_3097-1.txt jobOptions_KsKpi_data_3097-0.txt
           sed -i "s/ApplicationMgr\.EvtMax = -1/ApplicationMgr\.EvtMax = 20000/g" jobOptions_KsKpi_data_3097-0.txt
           sed -i "s/KsKpi_data_3097-1\.root/KsKpi_data_3097-0\.root/g" jobOptions_KsKpi_data_3097-0.txt
	       ;;

    0.1.3) echo "Test for data" 
           echo "have you changed test number?(yes / no)
           ./run/KsKpi/jobs_text/data/3097/jobOptions_KsKpi_data_3097-0.txt"
           read opt
           if [ $opt == "yes" ]
               then
               echo "now in yes"  
               cd ./run/KsKpi/jobs_text/data/3097
               boss.exe jobOptions_KsKpi_data_3097-0.txt
           else
               echo "Default value is 'no', please change test number."
           fi
           ;;

    0.1.4) echo "Submit Condor jobs on data ---- 2..." 
           cd ./run/KsKpi/gen_script
           ./sub_jobOption_file_Data_KsKpi.sh 705
	       ;;

    0.1.5) echo "Synthesize data root files..."
           cd ./run/KsKpi/gen_script
           ./syn_Data_KsKpi.sh 705
	       ;;

esac

}

sub_0_2() {

case $option in
    
    # --------------------------------------------------------------------------
    #  Signal MC  
    # --------------------------------------------------------------------------

    0.2.1) echo "Generate MC samples ---- Simulation && Reconstruction ..."
           echo "which MC sample do you want to simulate?"
           echo "    KsKpi       --> e+e- --> KsK+pi- + c.c."
           read opt
           if [ $opt == "KsKpi" ]; then
               cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/gen_script/gen_mc/KsKpi
               ./sub_jobOption_Ks_K_Pi.sh 705
           fi
	       ;;

    0.2.2) echo "Generate MC samples ---- Event Selection ..."
           echo "which MC sample do you want to select?"
           echo "    KsKpi       --> e+e- --> KsK+pi- + c.c."
           read opt
           if [ $opt == "KsKpi" ]; then
               cd /besfs/groups/cal/dedx/$USER/bes/KsKpi/run/KsKpi/gen_script/gen_mc/KsKpi
               ./sub_Ana_Ks_K_Pi.sh 705
           fi
	       ;;

    0.2.3) echo "Synthesize inclusive MC samples root files ..."
           cd ./run/KsKpi/gen_script
           echo "which MC sample do you want to synthesize?"
           echo "    KsKpi       --> e+e- --> KsK+pi- + c.c."
           read opt
           if [ $opt == "KsKpi" ]; then
               ./syn_MC_KsKpi.sh 705 KsKpi
           else
               echo "Please add the MC simulation joboption files!"
           fi
	       ;;

esac

}
case $option in
    
    # --------------------------------------------------------------------------
    #  Data  
    # --------------------------------------------------------------------------

    0.1) echo "Running on data sample..."
         usage_0_1 
         echo "Please enter your option: " 
         read option  
         sub_0_1 option 
	     ;;

    0.1.*) echo "Running on data sample..."
           sub_0_1 option  
           ;;  
        
    # --------------------------------------------------------------------------
    #  Signal MC 
    # --------------------------------------------------------------------------

    0.2) echo "Generating and running on signal MC sample..."
         usage_0_2
         echo "Please enter your option: " 
         read option  
         sub_0_2 option 
	     ;;

    0.2.*) echo "Generating and running on signal MC sample..."
           sub_0_2 option  
           ;;  

esac
