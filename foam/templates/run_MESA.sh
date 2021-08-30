#!/bin/bash

# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export Z_ini= $1
# export M_ini= $2
# export log_Dmix= $3
# export aov= $4
# export fov= $5
# export output_dir= $6
# export work_dir= $7

export WORK_DIR="${7}"
export LOG_DIR="${6}"
export OUTPUT_DIR=$LOG_DIR/Zini"$1"_Mini"$2"/

# Preparation:
mkdir -p $OUTPUT_DIR/history/
mkdir -p $OUTPUT_DIR/profiles/
mkdir -p $OUTPUT_DIR/gyre/
mkdir -p $OUTPUT_DIR/preMS/
cd $LOG_DIR


mkdir -p job_logs job_errs
log_file=$LOG_DIR/job_logs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5".log
err_file=$LOG_DIR/job_errs/Zini"$1"_Mini"$2"_logD"$3"_aov"$4"_fov"$5".err

cd $WORK_DIR

echo start: $(date)>>$log_file
# Execution:
./rn <<< $(echo \"$OUTPUT_DIR\", $1, $2, $3, $4, $5) 1>>$log_file  2>>$err_file

echo stop: $(date)>>$log_file

find $err_file -type f -size 0 -exec rm {} \;   # remove error file if empty
