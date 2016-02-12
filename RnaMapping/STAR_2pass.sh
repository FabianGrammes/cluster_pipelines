
test=1,0c66824bfbba50ee715658c4e1aeacf6fda7e7ff,1296,4234,194,1536,0
OIFS=$IFS;
IFS="|";


ta=($test)


# Join all SJ files
for ((i=0; i<${#SJa[@]}; ++i))
    do
	SJdir=${SJa[$i]}
	cat $SJdir/*SJ.out.tab >> star_2nd/SJ_all.tab
done

# Filter the joined file
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>1){print $1,$2,$3,strChar[$4]}}' star_2nd/SJ_all.tab > SJ_in.tab

echo '==>>FINISHED'
