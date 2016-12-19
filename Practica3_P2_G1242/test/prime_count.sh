#38106 seeds with the same result for 1024 bits
#128 iterations for the same seed with the same result for 1024 bits
#between 5 and 829 there are two (or one duplicate) which passes the test on one iteration
make
rm -f out.txt
for i in {944..1999}
do
	./primo -b 1024 -t 1 -s $i
	echo "$i"
	mv suc.txt suc$i.txt
	wc -l suc$i.txt | grep -P -o ^[1-9]+ >> out.txt
	rm suc$i.txt
done
