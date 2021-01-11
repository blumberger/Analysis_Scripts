fold="./100nsQuench/Clusters_Original"

cd $fold;

count=0
if [ -f "load_clusters.vmd" ] 
then
	rm load_clusters.vmd
fi

echo "display projection Orthographic" >> load_clusters.vmd
echo "axes location Off" >> load_clusters.vmd
for i in `find . -name "cluster*.xyz"`
do
	echo "mol new {$i} type {xyz} first 0 last -1 step 1 waitfor 1" >> load_clusters.vmd
	echo "mol modcolor 0 $count ColorID $count"         >> load_clusters.vmd     
	#echo "mol modstyle 0 $count Licorice 0.10000 27.000000 27.000000" >> load_clusters.vmd
	#echo "mol modmaterial 0 $count AOShiny" >> load_clusters.vmd
	#echo "mol modstyle 0 $count VDW 1.000000 12.000000" >> load_clusters.vmd

	count=$(( count + 1 ))
done

vmd -e load_clusters.vmd
