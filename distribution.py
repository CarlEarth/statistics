# The distribution includes the statistics module and example
import math # It provides some math function, such as exp(x), sin(x)...
import numpy as np
import matplotlib.pyplot as plt # Plot the figure
from scipy.integrate import quad # Integration
from scipy.optimize import fmin #Find the minimum

#********************************
## m-combination of a set S is a subset of n distinct elements of S.
## If the set has m elements, the number of n-combinations is
## equal to the binomial coefficient, c_mn
## (From Wiki: Combination)
def c_mn(m,n):
    a,b,c=1,1,1
    if (m < n or n < 0 or m < 0 ):
        print "fuction error"
        exit()
    else:
        if ((m-n) <= n):
            n=m-n
        for i in range(0, n):
            a=a*m
            b=b*c
            m-=1
            c+=1
    return a/b
#********************************
##solve the simple equation
##,which only has one simple cross point with f(x)=0 in the range (min,max).
## The algorithm is from Numerical analysis written by Richard L. Burden and J. Douglas Faires
def solve_eqn(f,min,max):
    i=1
    lk=min
    lj=max
    FA=f(lk)
    while(i < 10000):
        p = lk+(lj-lk)/2.0
        FP = f(p)
        if(FP == 0 or abs(FP) < 0.0000001):
            return p
        i=i+1
        if (FA*FP > 0.0):
            lk=p
            FA=FP
        else:
            lj=p
    print "error no~"
    return "SOLVE EQN ERROR"
#********************************
def Quad_f(x,M,S2): # Normal approximation to the log-likelihood: y=-(x-MLE)^2/(2*S^2) ~Xsquare
    a = - (x-M)**2.0/(2.0*S2) # Quadratic formula
    return a
#********************************
class binary() : #binomial distribution module
    def __init__(self, num_bin, success_p): #(number of bin, the probabolity of success)
        self.num = num_bin
        self.success_p = success_p
        self.E = num_bin*success_p #Expectation value
        self.var = num_bin*success_p*(1.0-success_p) #variance

    def func(self,num_success): # The probability of the "number of success"
        x=c_mn(self.num , num_success)*\
		(self.success_p**num_success)*\
		((1.-self.success_p)**(self.num-num_success))
        #P(n)=c_mn * p^n *(1-p)^(m-n)
        return x

    def interval(self,y,z):
        # The probability that the number of success events is between y and z. (y <= N < z)
        a,b=0.0,0.0
        for i in range(y, z):
            a=self.func(i)
            b=b+a
	
	return b
    def entropy(self,num_success): # The probability of the "number of success"
        x=c_mn(self.num , num_success)
        x=math.log(x)
        return x
	
#********************************

class like_binary() :#Likelihood of binomial distribution module
    
    def __init__(self, num_bin, num_success): # (number of bin, the number of success)
        self.num = num_bin
        self.num_success = num_success
        self.ext=solve_eqn(self.diff_ll,0.1,0.9)
	# find the maximum of the log-likelihood correponding to success_p
        # for calculting the log-likelihood ratio
    
    def l(self,success_p): # likelihood of the certain success probability
        x=c_mn(self.num , self.num_success)*\
            (success_p**self.num_success)*\
                ((1.-success_p)**(self.num-self.num_success))
        self.E = self.num *success_p
        self.var = self.num*success_p*(1.0-success_p)
        return x
    def ll(self,success_p): # log-likelihood of the certain success probability
        return math.log(self.l(success_p))

    def diff_ll(self,success_p): # derivative of log-likelihood with respect success probability
        h = 10**-10 #
        return (self.ll(success_p + h) - self.ll(success_p - h)) / (2 * h)
    
    # def extremum(self):
    #    a =solve_eqn(self.diff_ll,0.0001,1.0)
    #    return a

    def llr(self,success_p): #log-likelihood ratio
        return self.ll(success_p)-self.ll(self.ext)

    def interval(self,y): # To do some probability integral~
        aa,cc= quad(lambda a:self.l(a), 0.001, 1.0)
        c,d= quad(lambda a:self.l(a), self.ext - y, self.ext + y)
        return c/aa -0.95
    def confi95(self): #Calculate the confidence 95%
        bb= solve_eqn(self.interval,0.001,0.5)
        print  "the 95% interval is ",self.ext-bb,"~",self.ext+bb
        print  "the 95% interval is ",self.llr(self.ext-bb),"~",self.llr(self.ext+bb)
        return bb

#********************************
#********************************
class poisson() : 

    def __init__(self, mu):
	self.var=mu
	self.E=mu    
	self.sigma=self.var**0.5
    def func(self,x):
	b=self.E**x*math.exp(-self.E)
	c=1
	for i in range(0,x):		
		b=b/float(c)
		c+=1
	#b=(2.0*math.pi*self.var)**(-0.5)*math.exp(-(x-self.E)**2.0/(2.0*self.var))
        return b
    def interval(self,y,z):
	a,b=0.0,0.0
	for i in range(y, z):
		a=self.func(i)
		b=b+a
	
	return b
#********************************
class like_poisson() :

	def __init__(self, time, num_success, lambda_min, lambda_max): 

		self.num_success = num_success
		self.time = time
		self.min = lambda_min
		self.max = lambda_max
		self.ext = solve_eqn(self.diff_ll, lambda_min, lambda_max)
	def l(self ,ld):
		var =self.time *ld
		E =self.time * ld
		sigma = var**0.5
		b=E**self.num_success*math.exp(-E)
		c=1
		for i in range(0,self.num_success):
			b=b/float(c)
			c+=1
		#b=(2.0*math.pi*self.var)**(-0.5)*math.exp(-(x-self.E)**2.0/(2.0*self.var))
		return b
	def ll(self,ld):
		return math.log(self.l(ld))
    
	def diff_ll(self,ld):
		h = 0.0000001 # 0.0001
		return (self.ll(ld + h) - self.ll(ld - h)) / (2 * h)
    
	def llr(self,ld):
		a=self.ll(ld)
		b=self.ll(self.ext)
		return a-b

	def interval(self,y): # To do some probability integral~
        	aa,cc= quad(lambda a:self.l(a), self.min, self.max)
        	c,d= quad(lambda a:self.l(a), self.ext - y, self.ext + y)
        	return c/aa -0.95

	def confi95(self): #Calculate the confidence 95%
        	bb= solve_eqn(self.interval,0.0,(self.max-self.min)/2.0)
        	print  "the 95% interval is ",self.ext-bb,"~",self.ext+bb
        	print  "the 95% interval is ",self.llr(self.ext-bb),"~",self.llr(self.ext+bb)
        	return bb

#********************************
#********************************

class Gaussian() : 

    def __init__(self, mu, sigma2):
        self.var=sigma2
        self.E=mu
        self.sigma=self.var**0.5
    def func(self,x):
        b=(2.0*math.pi*self.var)**(-0.5)*math.exp(-(x-self.E)**2.0/(2.0*self.var))
        return b
    def interval(self,y,z):
	return quad(lambda a:self.func(a), y, z)

class like_Gaussian() : 

    def __init__(self, mu, sigma2):
        self.var=sigma2
        self.E=mu
        self.sigma=self.var**0.5
    def func(self,x):
        b=(2.0*math.pi*self.var)**(-0.5)*math.exp(-(x-self.E)**2.0/(2.0*self.var))
        return b
    def interval(self,y,z):
	return quad(lambda a:self.func(a), y, z)


