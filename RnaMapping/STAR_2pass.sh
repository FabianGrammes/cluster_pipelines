
test=1,0c66824bfbba50ee715658c4e1aeacf6fda7e7ff,1296,4234,194,1536,0
OIFS=$IFS;
IFS="|";


ta=($test)

for ((i=0; i<${#ta[@]}; ++i))
    do
echo ${ta[$i]}
done
