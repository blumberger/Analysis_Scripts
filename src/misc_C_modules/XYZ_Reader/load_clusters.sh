fold="./100nsQuench/Test"

cd $fold;

count=0
if [ -f "load_clusters.vmd" ] 
then
	rm load_clusters.vmd
fi
for i in `find . -name "cluster*.xyz"`
do
	echo "mol new {$i} type {xyz} first 0 last -1 step 1 waitfor 1" >> load_clusters.vmd
	echo "mol modcolor 0 $count ColorID $count"         >> load_clusters.vmd     
	#echo "mol modstyle 0 $count Licorice 0.150000 17.000000 17.000000" >> load_clusters.vmd
	#echo "mol modstyle 0 $count VDW 1.000000 12.000000" >> load_clusters.vmd

	count=$(( count + 1 ))
done

vmd -e load_clusters.vmd
