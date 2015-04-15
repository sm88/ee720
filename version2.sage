import random
import time
import multiprocessing

def indexCalculus(g,h,q,N,pid,dic):
    t0=time.time()
    
    #g generator
    #q modulus
    #h argument
    #N largest prime factor in factor base (r) factors  
   
    factor_base=[]
    
    i=2

    while(i<N):#N is the upper bound for the largest prime in factor baseT 
        #print i
        factor_base.append(i)
        i=next_prime(i)
    #print factor_base
    r=len(factor_base)
    
    N0=factor_base[-1] #N0 largest prime in factor base 

    countr=1;
    k=getQ(1,q-2)
    checkInd=False
    list_relation=[]
    list_relation2=[]
    
    # find r relations
    while(1):
        if(countr==r):
            break
        relation=[0]*(r+1)
        relation[r]=k
        
        z=sage.rings.integer.Integer(pow(g,k,q))
        f=factor(z)
        list_f=list(f)
        
        
        if len(list_f)==0 or list_f[-1][0]>N0:
            k=getQ(1,q-2)
            continue

        for prime_factor in list_f :
            idx=factor_base.index(prime_factor[0])
            relation[idx]=prime_factor[1]
        
        list_relation.append(relation)
        list_relation2.append(relation[:-1])

        if checkInd:
            matrix_relation=matrix(list_relation2)
            rnk=matrix_relation.rank()
            if rnk==countr+1:
                countr=countr+1
            else:
                list_relation.pop()
                list_relation2.pop()
       
        checkInd=True
        k=getQ(1,q-2)
    
    M=matrix(list_relation)
    #print M,'\n'
 	#reduced_M=M.echelon_form()
    #print reduced_M
    #print "rank =", M.rank()
    #print "no of relation=",M.nrows()
    
    log_factor_base=solveLSE(M,q-1)
    if len(log_factor_base) == 0:
        print 'time for iteration',time.time()-t0, N, pid
        dic[pid]=-1
        return
    #print "log_factor_base=",log_factor_base
    s=0
    fp=[0]*r
    while(1): 
        y=sage.rings.integer.Integer((pow(g,s,q)*(h%q))%q)
        f=factor(y)
        list_f=list(f)
        if len(list_f)>0 and list_f[-1][0]<=N0:
            for temp in list_f:
                idx=factor_base.index(temp[0])
                fp[idx]=temp[1]
            break
        s=s+1
    #print "fp=",fp
    
    
    i=0
    sum=0
    while(i<r):
        sum=(sum%(q-1)+((log_factor_base[i])*fp[i])%(q-1))%(q-1)
        i=i+1
    x=sum-s
    print 'time for iteration',time.time()-t0, N, pid
    dic[pid]=x%(q-1)
    #return x%(q-1)
    return

# function to solve linear system of equation modulus q
# input-> augmented matrix M 
# input-> q
def solveLSE(M,q):
    
    rref_M=M.echelon_form() #row reduced echelon form of the given matrix
    rows=rref_M.nrows()
    cols=rref_M.ncols()
    
    solution=[0]*rows
    i=1;

    while(i<=rows):
        r=rows-i
        c=cols-1-i

        j=1
        temp=rref_M[r][cols-1]
        sum=0
        while(j<i):
            sum=(sum%q+((rref_M[r][c+j])*solution[c+j])%q)%q
            j=j+1
        temp=(temp-sum)%q
        #print rref_M[r][c]
        #print "row,col",r,c
        val=rref_M[r][c]
        #print "val=",val
        #print "last line",rref_M[r]
        gcd1=gcd(val,q)
        if gcd1==1:

            solution[c]=(temp*(sage.rings.integer.Integer(pow(val,-1,q))))%q

        elif gcd1.divides(temp):
            a=sage.rings.integer.Integer(val/gcd1)
            b=sage.rings.integer.Integer(q/gcd1)
            #print "a=",a
            #print "b=",b
            temp=sage.rings.integer.Integer(temp/gcd1)
            #print "temp=",temp
            k=pow(a,-1,b)
            k=sage.rings.integer.Integer(k)
            #print "k=",k
            solution[c]=(temp*k)%q
            #print solution[c]
        else:
        	#NO SOLUTION
            #print "val=",val
            #print "gcd1=",gcd1
            #print "temp=",temp
            return []

        i=i+1

    #print solution

    return solution

def runIC(g,h,q,B0,Np,delB=50):
	print 
	print '          INPUT ', '\n'

	print '      prime :', q 
	print '  generator :', g
	print 'h=g^x mod q :', h, '\n'
	print '-------- SOLVING --------','\n'
	mgr=multiprocessing.Manager()
	dic=mgr.dict()
	jobs=[]
	isDone=False
	t0=time.time()
    
    #Np no of process
    #B0 initial values of the upper bound 
    #delb difference between two successive bounds
	B=range(B0, B0+Np*delB+1,delB)

	while not isDone:
		dic=mgr.dict()
		jobs=[]
		for i in range(Np):
			p=multiprocessing.Process(target=indexCalculus,args=(g,h,q,B[i],i,dic,))
			jobs.append(p)
			p.start()
			
		for p in jobs:
			p.join()
			
		for i in dic.keys():
			h1=sage.rings.integer.Integer(pow(g,dic[i],q))
			if h1==h:
				x=dic[i]
				bound=B[i]
				isDone=True
				break
				
		if not isDone:
			print 'inconsistent results, re-running'
	print '\n','--------- SOLVED ---------','\n'
	print '         x :',x
	print '         B :',bound
	print 'Total Time :', time.time()-t0    

def getQ(a,b,prime=False): # function for geting random number in the range [a,b]
    robj=random.SystemRandom(time.time())
    if prime:
        return next_prime(robj.randint(a,b))
    else:
        return robj.randint(a,b)

def genInput(n):# function for generating input

	# INPUT  : n               -> no. of digit
	# OUTPUT : list input[q,g] -> containing the prime q, and the multiplicative generator of the group g 
	
	input=[]
	rnd=random.SystemRandom(time.time())
	while True:
		q=next_prime(rnd.randint(1*(10**(n-1)),9*(10**(n-1))))
		if(q<(10^n)-1):
			break
	input.append(q)
	rings=Integers(q)
	g=rings.multiplicative_generator()
	input.append(g)

	return input
