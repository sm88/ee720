import random
import time
import multiprocessing
import os
import signal

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
        dic[pid]=(-2,time.time()-t0)
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
    dic[pid]=(x%(q-1),time.time()-t0)
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
	procs={}
    
    #Np no of process
    #B0 initial values of the upper bound 
    #delb difference between two successive bounds	
	B=range(B0, B0+Np*delB+1,delB)
	while not isDone:
		dic=mgr.dict()
		jobs=[]
		doCheck=True
		for i in range(Np):
			dic[i]=(-1,-1)
			p=multiprocessing.Process(target=indexCalculus,args=(g,h,q,B[i],i,dic,))
			procs[i]=p
			jobs.append(p)
			p.start()
			
		#poll for answer
		while doCheck:
			noGood=0
			for i in dic.keys():
				(ans,t) = dic[i]	#answer,time
				#is some answer found
				if ans != -1:
					#check validity
					h1=sage.rings.integer.Integer(pow(g,ans,q))
					if h1==h:
						x=i
						bound=B[i]
						isDone=True
						doCheck=False
						#kill other processes
						#print "killing"
						for pr in procs.keys():
							if pr!=i:
								while procs[i].is_alive():
									os.kill(procs[i].pid,signal.SIGTERM)
									time.sleep(0.1)
						break
					else:
						#how many wrong in poll
						noGood=noGood+1
					#all wrong
					if noGood == Np:
						doCheck=True
						break
				
		if not isDone:
			print 'inconsistent results, re-running'
	print '\n','--------- SOLVED ---------','\n'
	print 'ans :',dic[x][0]
	print 'B :',bound
	print 'best time:',dic[x][1]
	tt=time.time()-t0
	print 'Total Time :', tt

	return dic[x][1]



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

def sim(n1,n2,B0):
	timeList=[]
	with open('n_ver_t.data','a') as f:
		for i in range(n1,n2+1):
			x=getQ(10^4,10^6)		
			inp=genInput(i)
			q=inp[0]
			g=inp[1]
			h=sage.rings.integer.Integer(pow(g,x,q))
			t=runIC(g,h,q,B0,5)
			timeList.append(t)	
			f.write(str(i)+'\t'+str(t)+'\t'+str(g)+'\t'+str(h)+'\t'+str(q)+'\t'+str(x)+'\n')
			f.flush()
			os.fsync(f)
	print timeList
