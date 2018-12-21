import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
import matplotlib.animation as animation
import random

blue = (31/255., 119/255., 180/255.)
orange = (255/255., 127/255., 14/255.)
green = (44/255., 160/255., 44/255.)
red = (214/255., 39/255., 40/255.)

class Particles:

	pts = None 
	N = None
	num = None
	tf = None
	h = None
	n = 0 # current step

	bc = 10 # boundary condition goes from -bc to bc

	def __init__(self, N, num, tf):
		self.N = N
		self.num = num
		self.tf = tf
		self.pts = np.zeros((4,N,num)) # (x,y,vx,vy),t,particles
		self.h = (tf)/float(N)
		self.kin = np.zeros(N) # kinetic energy
		self.pot = np.zeros(N) # potential energy
		self.Etot = np.zeros(N) # total energy
		self.px = np.zeros(N)
		self.py = np.zeros(N)

	def setInitConditions(self):
		
		#self.pts[0][0][0] = 0.
		#self.pts[1][0][0] =	0.
		#self.pts[0][0][1] = 2**(1/6)
		#self.pts[1][0][1] = 0.

		# gridding
		x = 0
		y = 0

		for a in range(self.num):

			# Staggered particles
			#self.pts[0][0][a] = 2*a*(-1)**a    
			#self.pts[1][0][a] = 1.5*a*(-1)**a  
			#if self.pts[1][0][a] > self.bc:   
			#	self.pts[2][0][a] = 0.5 #(-1)**(a+1)*0.01 
			#else:
			#	self.pts[2][0][a] = -0.5
			#self.pts[3][0][a] = 0. #(-1)**a*0.01   
			
			# shear grid (100 particles)
			self.pts[0][0][a] = -9 + 2*x			
			self.pts[1][0][a] = 9 - 2*y
			self.pts[3][0][a] = 0.
			if y <= 4:
				self.pts[2][0][a] = 0.5
			else:
				self.pts[2][0][a] = 0.5
			self.px[0] += self.pts[2][0][a]
			self.py[0] += self.pts[3][0][a]
			x += 1
			if x > 9:
				x = 0
				y += 1

	def animate(self, i):

		for n in range(self.num):
			xn = np.zeros(int(self.N/step))
			yn = np.zeros(int(self.N/step))
			#for b in range(int(self.N/step)):
			#	xn[b] = x[step*b][n]
			#	yn[b] = y[step*b][n]
			xn[i-1] = x[i*step][n]
			yn[i-1] = y[i*step][n]
			line[n].set_xdata(xn[i-1])
			line[n].set_ydata(yn[i-1])

	def verlet(self):
		for i in range(self.num):
			Fx, Fy = self.Force(i ,self.n-1, pot = False)
			self.pts[2][self.n][i] = self.pts[2][self.n-1][i] + 0.5*self.h*Fx
			self.pts[3][self.n][i] = self.pts[3][self.n-1][i] + 0.5*self.h*Fy

		for i in range(self.num):
			self.pts[0][self.n][i] = self.pts[0][self.n-1][i] + self.h * self.pts[2][self.n][i]
			self.pts[1][self.n][i] = self.pts[1][self.n-1][i] + self.h * self.pts[3][self.n][i]
	
		for i in range(self.num):
			Fx, Fy = self.Force(i ,self.n, pot = True)
			self.pts[2][self.n][i] = self.pts[2][self.n][i] + 0.5*self.h*Fx
			self.pts[3][self.n][i] = self.pts[3][self.n][i] + 0.5*self.h*Fy
			
			# Boundary Conditions
			if self.pts[0][self.n][i] > self.bc:
				self.pts[0][self.n][i] -= 2*self.bc
			if self.pts[0][self.n][i] < -self.bc:
				self.pts[0][self.n][i] += 2*self.bc
			if self.pts[1][self.n][i] > self.bc:
				self.pts[1][self.n][i] -= 2*self.bc
			if self.pts[1][self.n][i] < -self.bc:
				self.pts[1][self.n][i] += 2*self.bc

	def Force(self, i, n, pot):
		Fx = 0
		Fy = 0 
		for a in range(num):
			if i != a:
				dx = (self.pts[0][n][a]-self.pts[0][n][i])
				dy = (self.pts[1][n][a]-self.pts[1][n][i])
				r = (self.pts[0][n][a]-self.pts[0][n][i])**2+(self.pts[1][n][a]-self.pts[1][n][i])**2

				# Nearest Neightbor Boundary Conditions
				if (dx > 0) and (dy > 0):
					if (dx > self.bc) and (dy <= self.bc):
						dx = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])
						r = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-self.pts[1][n][i])**2
					
					if (dx <= self.bc) and (dy > self.bc):
						dy = (self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-self.pts[0][n][i])**2+(self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])**2
					
					if (dx > self.bc) and (dy > self.bc):
						dx = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])
						dy = (self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])**2
				
				if (dx <= 0) and (dy <= 0):
					if (dx <= -self.bc) and (dy > -self.bc):
						dx = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])
						r = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-self.pts[1][n][i])**2
					
					if (dx >= -self.bc) and (dy <= -self.bc):
						dy = (self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-self.pts[0][n][i])**2+(self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])**2
					
					if (dx <= -self.bc) and (dy <= -self.bc):
						dx = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])
						dy = (self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])**2

				if (dx <= 0) and (dy > 0):
					if (dx <= -self.bc) and (dy <= self.bc):
						dx = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])
						r = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-self.pts[1][n][i])**2
					
					if (dx > -self.bc) and (dy > self.bc):
						dy = (self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-self.pts[0][n][i])**2+(self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])**2
					
					if (dx <= -self.bc) and (dy > self.bc):
						dx = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])
						dy = (self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]+2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-2*self.bc-self.pts[1][n][i])**2
				
				if (dx > 0) and (dy <= 0):
					if (dx > self.bc) and (dy > -self.bc):
						dx = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])
						r = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]-self.pts[1][n][i])**2
					
					if (dx <= self.bc) and (dy <= -self.bc):
						dy = (self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-self.pts[0][n][i])**2+(self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])**2
					
					if (dx > self.bc) and (dy <= -self.bc):
						dx = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])
						dy = (self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])
						r = (self.pts[0][n][a]-2*self.bc-self.pts[0][n][i])**2+(self.pts[1][n][a]+2*self.bc-self.pts[1][n][i])**2

				r = r**(1/2)
				if pot:
					self.pot[self.n] += 4*((1/r**12)-(1/r**6))
				Fx += -48*((1/r**14)-(1/2)*(1/r**8))*dx
				Fy += -48*((1/r**14)-(1/2)*(1/r**8))*dy
		return Fx, Fy
	
	def evolve(self):
		
		for i in range(1, self.N):
			self.n += 1
			self.verlet()
			for a in range(self.num):
				self.px[i] += self.pts[2][self.n][a]
				self.py[i] += self.pts[3][self.n][a]
			
		self.pot = self.pot/2
		for k in range(self.num):
			for n in range(self.N):
				self.kin[n] += 0.5*(self.pts[2][n][k]*self.pts[2][n][k]+self.pts[3][n][k]*self.pts[3][n][k])
		self.Etot = self.kin + self.pot

if __name__=='__main__':

	N = 5000 # max N
	num = 100 # number of particles
	tf = 50
	step = 5 # frame rate of animation

	p = Particles(N, num, tf)
	p.setInitConditions()
	p.evolve()
	#ax = p.plot()
	#p.plot(ax)

	time = np.linspace(0,tf,N)

	fig3, ax3 = plt.subplots( 1, 1)
	ax3.plot(time, p.kin, label='Kinetic')
	ax3.plot(time, p.pot, label='Potential')
	ax3.plot(time, p.Etot, label='Total')
	ax3.set_ylabel(r'$Energy$', fontsize=16)
	ax3.set_xlabel(r'$Time$', fontsize=16)
	ax3.set_title("Energy")
	ax3.legend( loc='best' )

	fig2, ax2 = plt.subplots( 1, 1)
	ax2.plot(time, p.px, label='Px')
	ax2.plot(time, p.py, label='Py')
	ax2.set_ylabel(r'$Momentum$', fontsize=16)
	ax2.set_xlabel(r'$Time$', fontsize=16)
	ax2.set_title("Momentum")
	ax2.legend( loc='best' )

	fig, ax = plt.subplots( 1, 1)
	t = np.linspace(1, 25, int(N/step))
	x = np.zeros((N,num))
	y = np.zeros((N,num))

	for n in range(num):
		for i in range(N):
			x[i][n] = p.pts[0][i][n]
			y[i][n] = p.pts[1][i][n]

	line = list()
	for n in range(num):
		line.append(ax.plot(x[:][n], y[:][n], lw=2, marker='o')[0])
	ax.set_xlim(-10, 10)
	ax.set_ylim(-10, 10)
	ax.set_xlabel(r'$X$', fontsize=16)
	ax.set_ylabel(r'$Y$', fontsize=16)

	anim = FuncAnimation(fig, p.animate, interval=1, frames=len(t)-1)
	#Writer = animation.writers['ffmpeg']
	#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
	#anim.save('/Users/nicholasfaucher/Documents/Physics/Computational_Physics/Final_Project/animation.mp4',writer=writer)

	plt.show()
