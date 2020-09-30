files_to_truncate="run-pseudo-hamilt-1.xyz run-pos-1.xyz run-coeff-1.xyz"
if [ -f "min_steps.sh" ]
then
   source min_steps.sh
fi

all_min_steps=""
for i in $files_to_truncate
do
   echo "Looking for $i files in run-fssh folders"

   var_name=$(echo `echo $i | grep '[a-z]*' -Poh` | tr " " "_")
   poss_min_step="${!var_name}"
   if [ "$poss_min_step" == "" ]
   then
      # Get the length of 1 step
      lines_in_step=`python -c "l='''$(grep "i *= *[0-9]*" run-fssh-0/$i -n | awk '{print $1}' | grep "[0-9]*" -Poh)'''.split(); print(int(l[1]) - int(l[0]))"`
   
      # First get num of steps
      all_steps=""
      all_files=`find run-fssh* -name "$i"`
      for fn in $all_files
      do
         all_steps="$all_steps `tail -$lines_in_step $fn | head -2 | grep "i *= *[0-9]*" -Poh | grep "[0-9]*" -Poh`"
      done
   
      # Get the minimum step in the list
      min_step=`python -c "print(min(map(int, '$all_steps'.split())))"`
      echo "Min Step: $min_step"
      
      echo "$var_name=$min_step" >> min_steps.sh

   else
      min_step=$poss_min_step
   fi

   all_min_steps="$all_min_steps $min_step"
done
min_step=`python -c "print(min(map(int, '$all_min_steps'.split())))"`

echo "Truncating files to the $min_step'th step"

for i in $files_to_truncate
do
   all_files=`find run-fssh* -name "$i"`
   lines_in_step=`python -c "l='''$(grep "i *= *[0-9]*" run-fssh-0/$i -n | awk '{print $1}' | grep "[0-9]*" -Poh)'''.split(); print(int(l[1]) - int(l[0]))"`

   # Get the line of the min step in each file and truncate
   for fn in $all_files
   do
      end_line_num=`grep "i = *$min_step" $fn -n | awk '{print $1}' | grep "[0-9]*" -Poh`
      echo "End Line: $end_line_num"
      end_line_num=`python -c "print($end_line_num-2 + $lines_in_step)"`
      start_line_num=`grep "= *0.000" $fn -n | tail -1 | awk '{print $1}' | grep "[0-9]*" -Poh`
      echo "Start Line: $start_line_num"
      start_line_num=`python -c "print($end_line_num - $start_line_num + 2)"`
      trunc_filename=`echo $fn | awk -F'.' '{print $1"_truncated.xyz"}'`

      echo "$trunc_filename:   Start Line: $start_line_num      End Line: $end_line_num"
      if [ -f $trunc_filename ]
      then
         rm $trunc_filename
      fi

      head -$end_line_num $fn | tail -$start_line_num > $trunc_filename
      echo "$trunc_filename complete (`wc -l $trunc_filename` lines)"
   done

done
