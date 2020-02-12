#!/bin/bash

# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export Z_ini= $1
# export M_ini= $2
# export log_Dmix= $3
# export aov= $4
# export fov= $5
# export Xc= $6
# GYRE inlist = $7
# export output_dir= $8

export OUTPUT_DIR="${8}"
export LOG_DIR="${OUTPUT_DIR%/*}"

# Preparation:
mkdir -p $OUTPUT_DIR
cd $LOG_DIR
mkdir -p job_errs job_logs

if [[ $7 == *"g_modes"* ]]; then
  log_file=$LOG_DIR/job_logs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6"_gmodes.log
  err_file=$LOG_DIR/job_errs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6"_gmodes.err
elif [[ $7 == *"p_modes"* ]]; then
  log_file=$LOG_DIR/job_logs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6"_pmodes.log
  err_file=$LOG_DIR/job_errs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6"_pmodes.err
else
  log_file=$LOG_DIR/job_logs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6".log
  err_file=$LOG_DIR/job_errs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5"_Xc"$6".err
fi

# Execution:
cd $OUTPUT_DIR
$GYRE_DIR/bin/gyre $7 1>>$log_file 2>>$err_file

find $err_file -type f -size 0 -exec rm {} \;   # remove error file if empty
