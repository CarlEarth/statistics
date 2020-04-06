import numpy as np
import matplotlib.pyplot as plt
#binomial distribution example
from distribution import binary, like_binary
from distribution import poisson, like_poisson
from distribution import Gaussian
import matplotlib.gridspec as gridspec
from math import log

#****************************
multipictue=0
Q1=0 #Q1: If X~Poi(30), please find P(x<=20). #Compare them with Gaussian 
Q2=0 #Q2: #one poisson eqn add another poisson eqn. Poi(30)+ Poi(30)=Poi(60)
     #Compare them with Gaussian 
Q3=1 #Additivity of entropy for two spin systems
#****************************

print "The following are the Poisson examples:"
#Show the Q1 result
if(Q1==1):
	tp1= poisson(30) #If X~Poi(30), please find P(x<=20).
	print "poisson distribution:"
	val1 = tp1.interval(0,21)
	print "Expectation value, variance=",tp1.E,tp1.var
	print "P(x<=20)=",val1
	#PAUSE
	raw_input("")
	tg3=Gaussian(tp1.E,tp1.var) 
	print "Use Gaussian to approximate this example:"
	val1,err = tg3.interval(0.0,20.5)
	print "E,var=",tg3.E,tg3.var
	print "P(x<=20.5)=",val1
	#PAUSE
	raw_input("")
	#plot
	prop1= []
	prog3= []
	k_range=range(0,60)
	for k in k_range:
		prop1.append(tp1.func(k))
		prog3.append(tg3.func(k))
	plt.figure()# If you need to plot a new figure, please add this line.
	plt.bar(k_range, prop1, color="g",
		     label="Poisson")
	plt.plot(k_range, prog3, color="b",
		     label="Gaussian")
	plt.legend(loc="best")
	plt.xlabel('k')
	plt.ylabel('probability')
	plt.show()

	


if(Q2==1):
	#one poisson eqn add another poisson eqn.
	# It make the the number of the events increase.
	# For the example, one gorup of particle Poi(30) add 
	#another group of particle Poi(30). Then, we'll get the P(60)
	tp1= poisson(30) #If X~Poi(30), please find P(x<=20).
	print "poisson distribution:"
	val1 = tp1.interval(0,41)
	print "Expectation value, variance=",tp1.E,tp1.var
	print "P(x<=40)=",val1
	#PAUSE
	raw_input("")
	tg3=Gaussian(tp1.E,tp1.var) 
	print "Use Gaussian to approximate this example:"
	val1,err = tg3.interval(0.0,40.5)
	print "E,var=",tg3.E,tg3.var
	print "P(x<=40.5)=",val1
	#PAUSE
	raw_input("")
	#plot
	prop1= []
	prog3= []
	k_range=range(0,120)
	for k in k_range:
		prop1.append(tp1.func(k))
		prog3.append(tg3.func(k))
	plt.figure()# If you need to plot a new figure, please add this line.
	plt.bar(k_range, prop1, color="g",
		     label="Poisson")
	plt.plot(k_range, prog3, color="b",
		     label="Gaussian")
	plt.legend(loc="best")
	plt.xlabel('k')
	plt.ylabel('probability')
	plt.show()
	
if(Q3==1):
	de2= like_poisson(160,8,0.00001,0.18)
	print "ext step"
	print "ext=",de2.ext
	de2.confi95()

	r_range=np.linspace(0.01,0.1,1000)
	loglikeratio=[]
	for l in r_range:
	    loglikeratio.append(de2.llr(l))

	plt.figure()
	plt.plot(r_range, loglikeratio, color="b",
		 label="poisson")
	plt.legend(loc="best")
	plt.xlabel('la,mda')
	plt.ylabel('likelihhood of total success:8')
	plt.savefig('poisson1')
	
