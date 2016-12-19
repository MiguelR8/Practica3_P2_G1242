make
rm -f out.txt
for i in {1000..15000000..1000}
do
	python "test/primeslessn.py" $i >> out.txt
done
