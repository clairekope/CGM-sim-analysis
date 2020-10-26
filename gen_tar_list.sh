if [ -f tar_list.txt ]; then
    rm tar_list.txt
fi

for dir in DD*/
do
    if [ ! -f ${dir:0:6}.tar.gz ]; then
	echo ${dir:0:6} >> tar_list.txt
    fi
done
