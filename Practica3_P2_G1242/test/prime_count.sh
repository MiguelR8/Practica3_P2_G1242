#38106 seeds with the same result for 1024 bits
#128 iterations for the same seed with the same result for 1024 bits
make
rm -f out.txt
for i in {1..10}
do
	./primo -b 1024 -t $i 
	echo "$i"
	mv suc.txt suc$i.txt
	wc -l suc$i.txt | grep -P -o ^[1-9]+ >> out.txt
	rm suc$i.txt
done
