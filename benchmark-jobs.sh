CFLAGS="-O3 -DNDEBUG -fopenmp"
if [ `uname -m` = amd64 ] || [ `uname -m` = x86_64 ]
then
	CFLAGS="$CFLAGS -mavx2"
fi

cc $CFLAGS -o benchmark *.c

rm -f benchmark-n*.out

for k in `seq 20`
do
	for i in `seq 1 15` `seq 16 2 62` `seq 64 4 128`
	do
		printf %02d/%03d\\n $k $i
		OMP_NUM_THREADS=$i ./benchmark -n$((1<<34)) >>benchmark-n$i.out
	done
done
