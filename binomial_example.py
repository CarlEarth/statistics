import numpy as np
import matplotlib.pyplot as plt
#binomial distribution example
from distribution import binary, like_binary
from distribution import Quad_f #Quadratic formula
import matplotlib.gridspec as gridspec
from math import log

#****************************
multipictue=0
Q1=0 #Q1: if X~Bin(100,0.5), find P(X>60)
Q2=0 #Q2: X~(bin=10, N of success=3) and Y~(bin=100, N of success=30)
     #Compare them in likelihood,log-likelihood, log-likelihood ratio
Q3=0 #Q3: Additivity of entropy for two spin systems
Q4=1 #Q4: Compare the log-likelihood ratio with the quadratic formula approximation:
     #    Example:Bin(100,success event:20) 
#****************************

print "The following are the binary examples:"

	#Show the Q1 result
if(Q1==1):
	#Q1: if X~Bin(100,0.5), find P(X>60)
	tb1=binary(100,0.5)

	print "##############Q1###################"
	print "Q1: if X~Bin(100,0.3)"
	print "Expectation value:",tb1.E
	print "Variance:",tb1.var
	#Sum P(X>60)
	print "P(X>60)=",tb1.interval(61,100)
	#PAUSE
	raw_input("Press Enter to continue...")
	#plot
	binary1= []
	k_range=range(30,70)

	for k in k_range:
		binary1.append(tb1.func(k))
	plt.figure()
	plt.bar(k_range, binary1, color="g")
	plt.plot(k_range, binary1, color="b")
	plt.legend(loc="best")
	plt.xlabel('k')
	plt.ylabel('probability')
#	plt.title('Binomial distribution')
	plt.show()


if(Q2==1):
	#Q2: X~(bin=10, N of success=3) and Y~(bin=100, N of success=30),
	#Compare them in likelihood,log-likelihood, log-likelihood ratio
	print "##############Q2###################"
	print "X~(bin=10, N of success=3) and Y~(bin=100, N of success=30)"
	de = like_binary(10,3)
	de1= like_binary(100,30)
	print "MLE in X=",de.ext
	print "MLE in Y=",de1.ext
	#PAUSE
	raw_input("Press Enter to continue...")
	r_range=np.linspace(0.2,0.4,100)
	r1_range=np.linspace(0.2,0.4,100)
	like= []
	loglike=[]
	loglikeratio=[]
	like1= []
	loglike1=[]
	loglikeratio1=[]
	l=0.0
	for l in r_range:
		like.append(de.l(l))
		like1.append(de1.l(l))
	for l in r1_range:
		loglike.append(de.ll(l))
		loglikeratio.append(de.llr(l))
		loglike1.append(de1.ll(l))
		loglikeratio1.append(de1.llr(l))
	if (multipictue==0):
		plt.figure()
		plt.plot(r_range, like, color="b",
			 label="Binary3/10")
		plt.plot(r_range, like1, color="g",
			 label="Binary30/100")
		plt.legend(loc="best")
		plt.xlabel('p of success')
		plt.ylabel('likelihhood')
		plt.title=('likelihhood')
		plt.show()
		
		plt.figure()

		plt.plot(r_range, loglike, color="b",
			 label="Binary3/10")
		plt.plot(r_range, loglike1, color="g",
			 label="Binary30/100")
		plt.legend(loc="best")
		plt.xlabel('p of success')
		plt.ylabel('log-likelihhood')
		plt.title=('loglikelihhood')
		plt.show()
		
		plt.figure()
		plt.plot(r_range, loglikeratio, color="b",
			 label="Binary3/10")
		plt.plot(r_range, loglikeratio1, color="g",
			 label="Binary30/100")
		plt.legend(loc="best")
		plt.xlabel('p of success')
		plt.ylabel('loglikelihhood ratio')
		plt.title=('loglikelihhood ratio')
		plt.show()


	if (multipictue==1):
		plt.figure()
		gs = gridspec.GridSpec(1, 3)
		ax1 = plt.subplot(gs[0, 0])
		ax1.plot(r_range, like, color="b",
			 label="Binary:3/10")
		ax1.plot(r_range, like1, color="g",
			 label="Binary:30/100")
		ax1.set_title('likelihood')

		ax2 = plt.subplot(gs[0, 1])
		ax2.plot(r1_range, loglike, color="b",
			 label="Binary:3/10")
		ax2.plot(r1_range, loglike1, color="g",
			 label="Binary:30/100")
		ax2.set_title('log-likelihood')

		ax3 = plt.subplot(gs[0, 2])
		ax3.plot(r1_range, loglikeratio, color="b",
			 label="Binary:3/10")
		ax3.plot(r1_range, loglikeratio1, color="g",
			 label="Binary:30/100")
		ax3.set_title('log-likelihoodratio')


		plt.show()
if(Q3==1):
	#Given two systems 
	N1=1000
	N2=1000
	spin1 = binary(N1,0.5)  
	spin2 = binary(N2,0.5)
	s_hat=0
	s=100
	total_entropy=0
	lss=0
	for s1 in range(0,s):
		lsdot=spin1.entropy(s1)+spin2.entropy(s-s1)
		if (total_entropy < lsdot):
			total_entropy = lsdot
			lss=s1
	print"max lss=",lss
	print"total entropy=", total_entropy

if(Q4==1):
	de = like_binary(100,20)
	#print de.E
	#print de.var
	de.confi95()

	#
	r_range=np.linspace(0.01,0.3,100)
	like= []
	loglike=[]
	loglikeratio=[]
	second=[]
	like1= []
	loglike1=[]
	loglikeratio1=[]
	l=0.0

	#Calculate the 2nd order equation to approximate the Binary llr.
	n=100.0
	n_test=20.0
	M=n_test/n
	p=M
	s2=p*(1-p)/n

	for l in r_range:
	    like.append(de.l(l))
	    loglike.append(de.ll(l))
	    loglikeratio.append(de.llr(l))
	    second.append(Quad_f(l,M,s2))
	#Plot likelihhood of Binary (100) with total success:20
	plt.figure()
	plt.title("likelihood in Binary20/100")
	plt.plot(r_range, like, color="b",
		 label="Binary20/100")
	plt.legend(loc="best")
	plt.xlabel('p of success')
	plt.ylabel('likelihhood of total success:20')
	plt.show()
	#Plot llr and quad_f approximation  and compare
	#We will find the approximation is good in some range.	
	plt.figure()
	plt.title("Log-likelihood ratio and quad_approx")
	plt.plot(r_range, loglikeratio, color="b",
		 label="Binary20/100")
	plt.plot(r_range, second, color="g",
		 label="approximation")
	plt.legend(loc="best")
	plt.xlabel('p of success')
	plt.ylabel('loglikelihhood ratio of total success:20')
	plt.show()


	
