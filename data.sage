import random
import time

rng=random.SystemRandom(time.time())
ls=[]
rings=[]
gens=[]
for i in range(0,5):
	ls.append(next_prime(rng.randint(3*(10**9),9*(10**9))))
	rings.append(Integers(ls[i]))
	gens.append(rings[i].multiplicative_generator())
	
	print ls[i], gens[i]
	
ls.append(next_prime(rng.randint(3*(10**19),9*(10**19))))
rings.append(Integers(ls[-1]))
gens.append(rings[-1].multiplicative_generator())

print ls[-1], gens[-1]
	
	
	
