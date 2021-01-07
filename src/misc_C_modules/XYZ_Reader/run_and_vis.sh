slice_num=`grep "Slice[0-9]" instruct.txt -Poh | grep "[0-9]" -Poh`

./run_slicer.sh
cp load_data.vmd 1nsQuenching/Slice$slice_num
cd 1nsQuenching/Slice$slice_num
#vmd -e load_data.vmd
~/Mercury/bin/mercury pos-init.xyz
