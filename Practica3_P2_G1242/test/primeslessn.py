import sys

if (len(sys.argv) < 2):
	print "Se necesita una cota superior"
	sys.exit()
try:
	top = int(sys.argv[1])
except ValueError:
	print "La cota debe ser un numero en base 10"
	sys.exit()

with open('test/primes1.txt') as f:
	count=0
	for l in f:
		for i in l.split(' '):
			try:
				if int(i) < top:
					count += 1
				else:
					print count
					sys.exit()
			except ValueError:
				pass
