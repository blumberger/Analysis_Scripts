"""
This script will take positions and velocities from the \$pos_mol_folder
and ammend CP2K input files to run MD simulations.

This will work with multiple threads (by starting a new process on a new thread
with bash i.e. cp2k -i run.inp &)

This version will give an estimated time of completion too.
"""



# Some parameters
nproc=8
min_step=0
save_data=false
pos_mol_folder="../../Pos_From_Previous"
CP2K_EXE="/scratch/mellis/flavoured-cptk/cp2k/exe/local/cp2k.sopt"
csv_headers="Step,Time,Kin,Temp,Pot,Etot,CPU,bonds,bends,urey-bradley,torsions,impropers,opbend,recip,real,self,neut,bonded,vdw"

#CP2K_EXE=cp2k
#CP2K_EXE="mpirun -n 3 /scratch/mellis/flavoured-cptk/cp2k/exe/local/cp2k.popt"

# Creat the csv file
echo "$csv_headers" > all_data.csv



get_ES_val() {
  # Get a number from a log file -for a given substring to search for.
  # The number must be on the same line as the search term, if multiple
  # lines fulfill the search criteria then multiple values will be returned.
  name=$1
  proc_num=$2

  grep "$name" RUN_$proc_num.log | grep "\-*[0-9][0-9]*\.*[0-9]*E*[-+]*[0-9]*" -Poh
}

escape_for_sed() {
  # Will print a string that has all special characters escaped with \
  KEYWORD="$1";
  printf '%s\n' "$KEYWORD" | sed -e 's/[]\/$*.^[]/\\&/g';
}
pos_mol_folder_escaped=`escape_for_sed "$pos_mol_folder"`

displaytime() {
  # Will convert seconds to hours, minutes, seconds
  # Will print a string %d days %H hours %M minutes %s seconds
  # If a non number is entered a python error will be returned
  local T=`python -c "print(int($1))"`
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}

# make the step data folder -to save data if requested
if [ "$save_data" == "true" ]
then
   data_fold="step_data"
   if [ -d "$data_fold" ] 
   then
      rm -r "$data_fold"
   fi  
   mkdir "$data_fold"
fi


# Get the step nums
all_steps=`python3 -c "import os, re; l=os.listdir('$pos_mol_folder'); print(' '.join(map(str, sorted(list(set([int(re.findall('[0-9    ]+', i)[0]) for i in l if 'py' not in i]))))))"`
declare -a step_nums
for i in $all_steps
do
   step_nums+=($i)
done

nsteps=${#step_nums[@]}
nstepsO=nsteps
nsteps=`python3 -c "print($nsteps - $min_step)"`

niter=`python3 -c "print(int($nsteps // $nproc) if $nsteps % $nproc == 0 else int($nsteps // $nproc) + 1)"`
niter=`python -c "print($niter - 1)"`


echo "Need to complete $niter iterations to finish $nsteps snapshots."
steps_left=$niter

# Ammend files for the parallel sims
for proc_num in `seq 0 $nproc`
do
   cp run.inp run_$proc_num.inp
   sed -i s/"PROJECT_NAME           run"/"PROJECT_NAME           run_$proc_num"/g  run_$proc_num.inp 
   sed -i s/"FORCE_EVAL.include"/"FORCE_EVAL_$proc_num.include"/g  run_$proc_num.inp 

   cp FORCE_EVAL.include FORCE_EVAL_$proc_num.include
   sed -i s/"TOPOLOGY.include"/"TOPOLOGY_$proc_num.include"/g FORCE_EVAL_$proc_num.include

   cp TOPOLOGY.include TOPOLOGY_$proc_num.include
done

# Run the sims
for iter in `seq 0 $niter`
do
   iter_procs=`python -c "print($nproc - 1)"`
   for proc_num in `seq 0 $iter_procs`
   do
      seq_num=`python -c "print(($iter * $nproc) + $proc_num + $min_step)"`
      step_num=${step_nums[$seq_num]}
      if ! [ -f "$pos_mol_folder/pos_$step_num.xyz" ]
      then
         continue
      fi

      echo $step_num
      
      sed s/"COORD_FILE_NAME.*"/"COORD_FILE_NAME   $pos_mol_folder_escaped\/pos_$step_num.xyz"/g TOPOLOGY_$proc_num.include -i
      sed s/"\@INCLUDE.*vel.*"/"\@INCLUDE   $pos_mol_folder_escaped\/vel_$step_num.include"/g FORCE_EVAL_$proc_num.include -i
      
      /usr/bin/time -o "time_$proc_num.txt" -f "%e" $CP2K_EXE run_$proc_num.inp &> RUN_$proc_num.log 2> run.err &

   done

   steps_left=`python -c "print($steps_left - 1)"`
   wait

   cat time_* >> all_times.txt
   avg_time=`awk '{x+=$1; next} END{print x/NR}' all_times.txt`
   time_left=`python -c "print($steps_left * $avg_time)"`
   echo "time left ($steps_left steps to go): $(displaytime `python -c "print(int($time_left))"`)"


   for proc_num in `seq 0 $iter_procs`
   do
      # Get step num
      seq_num=`python -c "print(($iter * $nproc) + $proc_num + $min_step)"`
      step_num=${step_nums[$seq_num]}

      # Get bonds, angle and urey bradley
      bond_line=`grep "BOND *= *" RUN_$proc_num.log`
      bond_ang_ub=`python3 -c "c=627.5094740631;s='''$bond_line'''.split(); print(f'{float(s[2])/c},{float(s[5])/c},{float(s[8])/c}')"`
      # Get torsions, impropers, opbends
      tors_line=`grep "TORSION *= *" RUN_$proc_num.log`
      tors_imp_ob=`python3 -c "c=627.5094740631;s='''$tors_line'''.split(); print(f'{float(s[2])/c},{float(s[5])/c},{float(s[8])/c}')"`
      echo $tors_imp_ob
      # Get Electrostatics
      recip=`get_ES_val "INITIAL GSPACE E" $proc_num`
      real=`get_ES_val "ES REAL ENERGY = " $proc_num`
      self=`get_ES_val "SELF ENERGY COR" $proc_num`
      neut=`get_ES_val "NEUT. BACKGROUND" $proc_num`
      bonded=`get_ES_val "BONDED CORRECTION" $proc_num`
      vdw=`get_ES_val "LJ ENERGY = " $proc_num`


      # Add data to csv
      line=`tail -1 run_$proc_num-1.ener`
      new_line=`python -c "words = '$line'.split(); words[0] = '$step_num'; print(','.join(words))"`
      new_line="$new_line,$bond_ang_ub,$tors_imp_ob,$recip,$real,$self,$neut,$bonded,$vdw"
      echo $new_line >> all_data.csv
      
      if [ "$save_data" == "true" ]
      then
         mv run_$proc_num-pos-1.xyz "$data_fold/pos_$step_num.xyz"
         echo "run_$proc_num-vel-1.xyz,$step_num" | python3 trans_vel.py
         rm "run_$proc_num-vel-1.xyz"
      fi 
   done
done


# Clear up
for proc_num in `seq 0 $nproc`
do
   rm -f run_$proc_num.inp TOPOLOGY_$proc_num.include FORCE_EVAL_$proc_num.include run_$proc_num-pos-1.xyz run_$proc_num-vel-1.xyz run_$proc_num-1.ener
done

