#!/bin/sh
export SRC_DIR=/clusterfs/esd/esd2/gbisht/Projects/ACME/repo/ACME
export CLM_USRDAT_DOMAIN=domain_82x1_sparse_grid_c150528.nc
export CLM_USRDAT_SURDAT=surfdata_82x1_sparse_grid_c150528.nc
export CLM_USRDAT_DIR=/clusterfs/esd/esd2/gbisht/Projects/ACME/repo/misc-scripts-for-acme/matlab-script-for-sparse-grid/82x1_sparse_grid

cd $SRC_DIR/scripts

export GIT_HASH=`git log -n 1 --pretty=%h`
export RES=f19_g16
export RES=CLM_USRDAT
export COMPSET=I1850CRUCLM45BGC
export MAC=lawrencium-lr3
export COMPILER=intel
export PROJ=acme
export CASE_NAME=${RES}.${COMPSET}.${MAC}.${COMPILER}.${GIT_HASH}.`date +"%Y-%m-%d"`

./create_newcase -case ${CASE_NAME} -res ${RES} -compset ${COMPSET} -mach ${MAC} -project ${PROJ} -compiler ${COMPILER}


cd ${CASE_NAME}

./xmlchange -file env_run.xml -id DIN_LOC_ROOT -val /clusterfs/esd/esd2/gbisht/ccsm_inputdata
./xmlchange -file env_run.xml -id DIN_LOC_ROOT_CLMFORC -val /clusterfs/esd/esd2/gbisht/ccsm_inputdata/atm/datm7
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val 1901

./xmlchange -file env_run.xml -id ATM_DOMAIN_FILE -val ${CLM_USRDAT_DOMAIN}
./xmlchange -file env_run.xml -id ATM_DOMAIN_PATH -val ${CLM_USRDAT_DIR}
./xmlchange -file env_run.xml -id LND_DOMAIN_FILE -val ${CLM_USRDAT_DOMAIN}
./xmlchange -file env_run.xml -id LND_DOMAIN_PATH -val ${CLM_USRDAT_DIR}

cat >> user_nl_clm << EOF
fsurdat = '${CLM_USRDAT_DIR}/${CLM_USRDAT_SURDAT}'
EOF

./cesm_setup

# Check is using intel/2015.0.090
intel_2015_0_090=`cat env_mach_specific  | grep 'module load intel/2015.0.090'`
if [ -z "$intel_2015_0_090" ]; then
  perl -w -i -p -e 's@module load intel@module load intel/2015.0.090@' env_mach_specific
fi

perl -w -i -p -e 's@#SBATCH --account=acme@#SBATCH --account=lr_esd2@' ${CASE_NAME}.run
perl -w -i -p -e 's@#SBATCH --qos=lr_normal@#SBATCH --qos=condo_esd2@' ${CASE_NAME}.run

./${CASE_NAME}.build

#sbatch ${CASE_NAME}.run

