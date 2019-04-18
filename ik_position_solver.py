"""
Closed form solutions to solve the Inverse Kinematics of Delta Robot

INPUT: Goal position
OUTPUT: Desired thetas
"""

from math import sqrt
import numpy as np
pi = np.pi
arctan = np.arctan
sin = np.sin
cos = np.cos

RAD2DEG = 180.0/pi

class position(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

class deltaSolver(object):
	# Using dimensions from the Delta Robot paper examples (pg. 16)
	def __init__(self, sb = 0.567, sp = 0.076, L = 0.524, l = 1.244, h = 0.131, tht0 = (0, 0, 0)):
		(self.currTheta1, self.currTheta2, self.currTheta3) = tht0
		self.vel1 = 0
		self.vel2 = 0
		self.vel3 = 0

		self.sb = sb
		self.sp = sp
		self.L = L
		self.l = l
		self.h = h

		self.wb = (sqrt(3)/6) * self.sb
		self.ub = (sqrt(3)/3) * self.sb
		self.wp = (sqrt(3)/6) * self.sp
		self.up = (sqrt(3)/3) * self.sp
		
		self.a = self.wb - self.up
		self.b = self.sp/2 - (sqrt(3)/2) * self.wb
		self.c = self.wp - self.wb/2

	def solveTheta1(self, position):
		#Takes in an argument that is position class
		#Solves for Theta1
		E1 = 2*self.L*(position.y+self.a)
		F1 = 2*position.z*self.L
		G1 = position.x**2 + position.y**2 + position.z**2 + self.a**2 + self.L**2 + 2*position.y*self.a - self.l**2

		return self.angleSolver(E1, F1, G1, 1)

	def solveTheta2(self, position):
		E2 = -self.L*(sqrt(3) * (position.x + self.b) + position.y + self.c)
		F2 = 2*position.z*self.L
		G2 = position.x**2 + position.y**2 + position.z**2 + self.b**2 + self.c**2 + self.L**2 + 2*(position.x*self.b + position.y*self.c) - self.l**2

		return self.angleSolver(E2, F2, G2, 2)

	def solveTheta3(self, position):
		E3 = self.L * (sqrt(3) * (position.x - self.b) - position.y - self.c)
		F3 = 2*position.z*self.L
		G3 = position.x**2 + position.y**2 + position.z**2 + self.b**2 + self.c**2 + self.L**2 + 2*(-position.x*self.b + position.y*self.c) - self.l**2

		return self.angleSolver(E3, F3, G3, 3)

	def angleSolver(self, E, F, G, thetaID):
		t1 = (-F + sqrt(E**2 + F**2 - G**2))/(G - E)
		t2 = (-F - sqrt(E**2 + F**2 - G**2))/(G - E)
		thetaPossible1 = 2*arctan(t1)
		thetaPossible2 = 2*arctan(t2)

		if(thetaID == 1):
			currTheta = self.currTheta1
		
		elif(thetaID == 2):
			currTheta = self.currTheta2

		elif(thetaID == 3):
			currTheta = self.currTheta3

		# Calculate the difference between the possible angles that solves the quadratic with current angle. 
		thetaDiff1 = thetaPossible1 - self.currTheta1
		thetaDiff2 = thetaPossible2 - self.currTheta2

		# Return the theta that is closest to the current theta
		if(abs(thetaDiff1) < abs(thetaDiff2)):
			return thetaPossible1 * RAD2DEG
		else:
			return thetaPossible2 * RAD2DEG

	def ik(self, goal):
		return [self.solveTheta1(goal), self.solveTheta2(goal), self.solveTheta3(goal)]


# TEST CODE BELOW
ds = deltaSolver()
test = position(0, 0, -0.9)
print(test.x, test.y, test.z)
print(ds.ik(test))

test2 = position(0.3, 0.5, -1.1)
print(test2.x, test2.y, test2.z)
print(ds.ik(test2))