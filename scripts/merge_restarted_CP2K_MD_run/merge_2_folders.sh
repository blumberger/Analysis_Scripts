orig_fold="."
fold_to_merge="restart_609000"
out_fold="merged_data"


if ! [ -d "$out_fold" ]
then
   mkdir $out_fold
else
   rm -rf $out_fold/*
fi


# Do the pos and vel files
for typ in "pos" "vel"
do
   file1="$orig_fold/run-$typ-1.xyz"
   file2="$fold_to_merge/run-$typ-1.xyz"

   num_at=`head -1 $file1`
   num_lines_in_step=`python -c "print($num_at + 2)"`
   last_time_str=`tail -$num_lines_in_step $file1 | head -2 | grep -Poh "time *= *[0-9]*\.[0-9]*"`
   last_step_str=`tail -$num_lines_in_step $file1 | head -2 | grep -Poh "i *= *[0-9]*"`
   last_time=`echo $last_time_str  | awk '{print $3}'`
   last_step=`echo $last_step_str  | awk '{print $3}'`
   
   all_time_strings=`grep -Poh "time *= *[0-9]*\.[0-9]*" $file2 | awk '{print $3}'`
   new_time_strings=`python -c "print(' '.join(map(str, [float(i) + $last_time for i in '''$all_time_strings'''.split()])))"`
   
   # Get num of nums that need to change
   LEN_LIST=0
   for i in $all_time_strings
   do
      LEN_LIST=`python -c "print($LEN_LIST + 1)"`
   done
   
   # Remove the first step of the new data (as it is repeated)
   num_lines_f2=`wc -l $file2 | awk '{print $1}'`
   l=`python -c "print($num_lines_f2 - $num_lines_in_step)"`
   tail -$l $file2 > tmp.test
   num_lines_tmp=`wc -l tmp.test | awk '{print $1}'`

   # Loop over all items
   for i in `seq 0 $(python -c "print($LEN_LIST - 2)")`
   do
      old_time_str=`python -c "print('''$all_time_strings'''.split()[$i+1])"`
      new_time_str=`python -c "print($old_time_str + $last_time)"`
   
      line_to_change=`python -c "print(($i * $num_lines_in_step) + 2)"`
      line_1m_to_change=`python -c "print($line_to_change - 1)"`
      tail_start=`python -c "print($num_lines_tmp - $line_to_change)"`

      # Handle the timestep
      keyword=$(printf '%s\n' "$old_time_str" | sed -e 's/[]\/$*.^[]/\\&/g')
      escaped_replace=$(printf '%s\n' "$new_time_str" | sed -e 's/[]\/$*.^[]/\\&/g')
      replace_line=`head -$line_to_change tmp.test | tail -1`
      replace_line=`echo $replace_line | sed s/"$keyword"/"$escaped_replace"/`
      step_num=`echo $replace_line | grep -Poh "i *= *[0-9]*" | awk '{print $3}'`
      new_step_num=`python -c "print($step_num + $last_step)"`
      replace_line=`echo $replace_line | sed "s/i *= *[0-9]*/i = $new_step_num/"`

      head -$line_1m_to_change tmp.test > new_file.tmp
      echo $replace_line >> new_file.tmp
      tail -$tail_start tmp.test >> new_file.tmp

      mv new_file.tmp tmp.test

   
      printf "\rCompleted Step $i          \r"
   done
   
   cat $file1 tmp.test > $out_fold/run-$typ-1.xyz
   rm tmp.test

   echo "Done $typ                             "
done


# The energy file
python3 merge_ener.py $fold_to_merge $orig_fold
cat run-1.ener new_file.tmp > $out_fold/run-1.ener
echo "Done ener"



echo "Finished"
