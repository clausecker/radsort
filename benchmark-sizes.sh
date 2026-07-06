CFLAGS="-O3 -DNDEBUG -fopenmp"
if [ `uname -m` = amd64 ] || [ `uname -m` = x86_64 ]
then
	CFLAGS="$CFLAGS -mavx2"
fi

cc $CFLAGS -o benchmark *.c

rm -f benchmark-*.out

export OMP_NUM_THREADS=4

for k in `seq 20`
do
	for i in `seq 1 30`
	do
		printf %02d/%d\\n $k $((1<<i))
		./benchmark -n $((1<<i)) >>benchmark-$((1<<i)).out
		printf %02d/%d\\n $k $(((1<<i-1)*3))
		./benchmark -n $(((1<<i-1)*3)) >>benchmark-$(((1<<i-1)*3)).out
	done
done

