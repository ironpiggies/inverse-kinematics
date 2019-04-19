"""
Closed form solutions to solve the Inverse Kinematics of Delta Robot

INPUT: Goal position
OUTPUT: Desired thetas
"""

from math import sqrt
from scipy.optimize import fsolve
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time

pi = np.pi
arctan = np.arctan
sin = np.sin
cos = np.cos

RAD2DEG = 180.0/pi
DEG2RAD = pi/180.0

class position(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z


def rotz(theta):
	return np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta),  np.cos(theta), 0], [0, 0, 1]])

def roty(theta):
	return np.array([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])

def rotx(theta):
	return np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta),  np.cos(theta)]])


class deltaSolver(object):
	# Using dimensions from the Delta Robot paper examples (pg. 16)
	#def __init__(self, sb = 567, sp = 76, L = 524, l = 1244, h = 131, tht0 = (0, 0, 0)):
	def __init__(self, sb = 2*109.9852, sp = 109.9852, L = 304.8, l = 609.5144, h = 42.8475, tht0 = (0, 0, 0), swivel_limit = 10*DEG2RAD):
	# def __init__(self, sb = 0.567, sp = 0.076, L = 0.524, l = 1.244, h = 0.131, tht0 = (0, 0, 0)):
		self.swivel_limit = swivel_limit
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

		(xx, yy, zz)=self.FK((self.currTheta1, self.currTheta2, self.currTheta3))
		self.x = xx
		self.y = yy
		self.z = zz
		self.endpt = (self.x, self.y, self.z)
		(th1, th2, th3) = self.ik(position(self.x, self.y, self.z))
		self.thts = (th1, th2, th3)
		self.fig = plt.figure()
		# (xx, yy, zz)=self.FK((self.currTheta1, self.currTheta2, self.currTheta3))
		# self.x = xx
		# self.y = yy
		# self.z = zz
		# self.endpt = (self.x, self.y, self.z)
		# self.pos = position(self.x, self.y, self.z)
		# (th1, th2, th3) = self.ik(self.pos)
		# self.thts = (th1, th2, th3)
		# self.fig = plt.figure()

		self.plot(position(xx,yy,zz))

	def solveTheta1(self, position):
		#Takes in an argument that is position class
		#Solves for Theta1
		# print('z2', position.z**2)
		# print('part of g', self.a**2 + self.L**2 + 2*position.y*self.a - self.l**2)
		E1 = 2*self.L*(position.y+self.a)
		F1 = 2*position.z*self.L
		G1 = position.x**2 + position.y**2 + position.z**2 + self.a**2 + self.L**2 + 2*position.y*self.a - self.l**2

		# Singularity condition because leg is aligned with y-axis
		if(G1 - E1 == 0):
			return 0
		# print('E1', E1, 'F1', F1, 'G1', G1)
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
		# print(G)
		# print(E)
		# print(thetaID)
		if(thetaID == 1):
			currTheta = self.currTheta1
		
		elif(thetaID == 2):
			currTheta = self.currTheta2

		elif(thetaID == 3):
			currTheta = self.currTheta3

		# Calculate the difference between the possible angles that solves the quadratic with current angle. 
		thetaDiff1 = thetaPossible1 - currTheta
		thetaDiff2 = thetaPossible2 - currTheta
		# Return the theta that is closest to the current theta
		if(abs(thetaDiff1) < abs(thetaDiff2)):
			return thetaPossible1
		else:
			return thetaPossible2

	def ik(self, goal):
		return [self.solveTheta1(goal), self.solveTheta2(goal), self.solveTheta3(goal)]

	def FK(self,thts):
		# Works regardless of length unit. Angle units are in radians. 
		th1, th2, th3 = thts
		def simulEqns(inp):
			(x, y, z) = inp
			l = self.l
			L = self.L
			a = self.a
			b = self.b
			c = self.c
			eq1 = 2*z*L*sin(th1) + x*x + y*y + z*z - l*l + L*L + a*a + 2*y*a + 2*L*(y+a)*cos(th1) 
			eq2 = 2*z*L*sin(th2) + x*x + y*y + z*z - l*l + L*L + b*b + c*c + 2*x*b + 2*y*c - L*(sqrt(3)*(x+b)+y+c)*cos(th2) 
			eq3 = 2*z*L*sin(th3) + x*x + y*y + z*z - l*l + L*L + b*b + c*c - 2*x*b + 2*y*c + L*(sqrt(3)*(x-b)-y-c)*cos(th3)
			return (eq1, eq2, eq3)
		return fsolve(simulEqns,(0,0,-500))

	def vel(self, inputVec):
		#inputVec is a vector in meters per second. [xdot, ydot, zdot]
		
		a11 = self.x
		a12 = self.y+self.a + self.L*cos(self.currTheta1)
		a13 = self.z + self.L*sin(self.currTheta1)

		a21 = 2*(self.x+self.b) - sqrt(3)*self.L*cos(self.currTheta2)
		a22 = 2*(self.y + self.c) - self.L*cos(self.currTheta2)
		a23 = 2*(self.z + self.L*sin(self.currTheta2))

		a31 = 2*(self.x - self.b) + sqrt(3)*self.L*cos(self.currTheta3)
		a32 = 2*(self.y + self.c) - self.L*cos(self.currTheta3)
		a33 = 2*(self.z + self.L*sin(self.currTheta3))

		A = np.matrix([[a11, a12, a13], 
					   [a21, a22, a23], 
					   [a31, a32, a33]])

		b11 = self.L*((self.y + self.a)*sin(self.currTheta1) - self.z*cos(self.currTheta1))
		b22 = -self.L*((sqrt(3)*(self.x + self.b) + self.y + self.c)*sin(self.currTheta2) + 2*self.z*cos(self.currTheta2))
		b33 = self.L*((sqrt(3)*(self.x-self.b) - self.y - self.c)*sin(self.currTheta3) - 2*self.z*cos(self.currTheta3))
		
		B = np.matrix([[b11, 0, 0], 
					   [0, b22, 0], 
					   [0, 0, b33]])

		AB = np.matmul(A,B)
		thetadot = np.matmul(AB, inputVec)

		#return a vector of the angle velocities. [omega1, omega2, omega3]
		return thetadot

	def plot(self, pos = position(0, 0, -500)):
		(x, y, z) = (pos.x, pos.y, pos.z)
		thts = self.ik(pos)

		ax = self.fig.add_subplot(111, projection='3d')
		self.ax = ax
		ax.set_xlim3d(-400, 400)
		ax.set_ylim3d(-400, 400)
		ax.set_zlim3d(-900, 100)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.ion()
		plt.show()
		ax.set_xlabel('X [mm]')
		ax.set_ylabel('Y [mm]')
		ax.set_zlabel('Z [mm]')
		#Draw Origin
		ax.scatter(0,0,0, marker = '+', c = 'k')

		#Draw Base
		base1 = np.matrix([-self.sb/2,-self.wb,0]).transpose()
		base2 = np.matrix([self.sb/2,-self.wb,0]).transpose()
		base3 = np.matrix([0,self.ub,0]).transpose()
		basePts = np.hstack((base1, base2, base3, base1))
		basePts = np.array(basePts)
		ax.plot(basePts[0,:] ,basePts[1,:], basePts[2,:],c='k')
		
		#Plot Origin Axes
		p = np.array([0, 0, 0])
		a1 = p+np.array([100,0,0])
		a2 = p+np.array([0,100,0])
		a3 = p+np.array([0,0,100])
		ax.plot([p[0], a1[0]], [p[1], a1[1]], [p[2], a1[2]], c='r', marker = '<')
		ax.plot([p[0], a2[0]], [p[1], a2[1]], [p[2], a2[2]], c='g', marker = '<')
		ax.plot([p[0], a3[0]], [p[1], a3[1]], [p[2], a3[2]], c='b', marker = '<')

		#Plot End Platform
		p = np.array([[x, y, z]]).T
		BTp1 = p+np.array([[0, -self.up, 0]]).T
		BTp2 = p+np.array([[self.sp/2, self.up, 0]]).T
		BTp3 = p+np.array([[-self.sp/2, self.up, 0]]).T
		BTp = np.array(np.hstack((BTp1, BTp2, BTp3, BTp1)))
		self.myPts = ax.plot(BTp[0,:], BTp[1,:], BTp[2,:],c='darkviolet')

		#Plot linkages
		pt1B = np.array([[0,-self.wb,0]]).T
		pt1J = pt1B+np.array([[0, -self.L*cos(-thts[0]), self.L*sin(-thts[0])]]).T
		pt1P = BTp1
		pt2B = np.dot(rotz(2*np.pi/3), pt1B)
		pt2J = pt2B+np.dot(rotz(2*np.pi/3), np.array([[0, -self.L*cos(-thts[1]), self.L*sin(-thts[1])]]).T)
		pt2P = BTp2
		pt3B = np.dot(rotz(4*np.pi/3) , pt1B)
		pt3J = pt3B+np.dot(rotz(4*np.pi/3), np.array([[0, -self.L*cos(-thts[2]), self.L*sin(-self.thts[2])]]).T)
		pt3P = BTp3
		self.link1 = ax.plot([pt1B[0][0], pt1J[0][0], pt1P[0][0]], [pt1B[1][0], pt1J[1][0], pt1P[1][0]], [pt1B[2][0], pt1J[2][0], pt1P[2][0]], c='dimgrey')
		self.link2 = ax.plot([pt2B[0][0], pt2J[0][0], pt2P[0][0]], [pt2B[1][0], pt2J[1][0], pt2P[1][0]], [pt2B[2][0], pt2J[2][0], pt2P[2][0]], c='dimgrey')
		self.link3 = ax.plot([pt3B[0][0], pt3J[0][0], pt3P[0][0]], [pt3B[1][0], pt3J[1][0], pt3P[1][0]], [pt3B[2][0], pt3J[2][0], pt3P[2][0]], c='dimgrey')

		#Update the Figure
		self.fig.canvas.draw_idle()
		plt.pause(0.0001)

	def updatePlot(self, pos = position(0, 0, -500), color = None):
		if(color == None):
			(x, y, z) = (pos.x, pos.y, pos.z)
			thts = self.ik(pos)

			#Plot End Platform
			p = np.array([[x, y, z]]).T
			BTp1 = p+np.array([[0, -self.up, 0]]).T
			BTp2 = p+np.array([[self.sp/2, self.up, 0]]).T
			BTp3 = p+np.array([[-self.sp/2, self.up, 0]]).T
			BTp = np.array(np.hstack((BTp1, BTp2, BTp3, BTp1)))
			self.updateThings(self.myPts, BTp[0,:], BTp[1,:], BTp[2,:])

			#Plot linkages
			pt1B = np.array([[0,-self.wb,0]]).T
			pt1J = pt1B+np.array([[0, -self.L*cos(-thts[0]), self.L*sin(-thts[0])]]).T
			pt1P = BTp1
			pt2B = np.dot(rotz(2*np.pi/3), pt1B)
			pt2J = pt2B+np.dot(rotz(2*np.pi/3), np.array([[0, -self.L*cos(-thts[1]), self.L*sin(-thts[1])]]).T)
			pt2P = BTp2
			pt3B = np.dot(rotz(4*np.pi/3) , pt1B)
			pt3J = pt3B+np.dot(rotz(4*np.pi/3), np.array([[0, -self.L*cos(-thts[2]), self.L*sin(-thts[2])]]).T)
			pt3P = BTp3
			self.updateThings(self.link1, [pt1B[0][0], pt1J[0][0], pt1P[0][0]], [pt1B[1][0], pt1J[1][0], pt1P[1][0]], [pt1B[2][0], pt1J[2][0], pt1P[2][0]])
			self.updateThings(self.link2, [pt2B[0][0], pt2J[0][0], pt2P[0][0]], [pt2B[1][0], pt2J[1][0], pt2P[1][0]], [pt2B[2][0], pt2J[2][0], pt2P[2][0]])
			self.updateThings(self.link3, [pt3B[0][0], pt3J[0][0], pt3P[0][0]], [pt3B[1][0], pt3J[1][0], pt3P[1][0]], [pt3B[2][0], pt3J[2][0], pt3P[2][0]])

			#Update the Figure
			self.fig.canvas.draw_idle()
			plt.pause(0.1)
		else:
			self.ax.scatter(xs = pos[0], ys = pos[1],zs = pos[2], c = color)

	def update_lines(self, num, dataLines, lines) :
		for line, data in zip(lines, dataLines) :
			# note: there is no .set_data() for 3 dim data...
			line.set_data(data[0:2, num:num+2])
			line.set_3d_properties(data[2,num:num+2])
		return lines

	def updateThings(self, linesObj, xPts, yPts, zPts):
		linesObj[0].set_data(xPts, yPts)
		linesObj[0].set_3d_properties(zPts)

	def check_constraints(self, motorID, endpos, theta1):
		#endpos is a tuple

		#assign lengths
		length_to_motor = self.wb
		length_to_attach = self.up

		#create an array of the end position
		xe = np.array([endpos[0], endpos[1], endpos[2]])

		#all rotation matrices along the z axis
		if(motorID == 1):
			#rotation matrix is rotation -90 degrees
			rotation_matrix = rotz(-pi/2)
		if(motorID == 2):
			#rotation_matrix is rotation 120-90 degrees
			rotation_matrix = rotz(2*pi/3-pi/2)
		if(motorID == 3):
			#rotation matrix is rotation 240-90 degrees
			rotation_matrix = rotz(4*pi/3-pi/2)

		#rotate the end effector position so that the math is easy
		xe = np.dot(rotation_matrix,xe)

		#Find the end position of the attachment
		bvec = xe - np.array([0, length_to_attach, 0])
		x0 = np.array([length_to_motor, 0, 0])
		#Solve two equations to solve for theta3
		theta2 = np.arcsin((bvec[2] - x0[2] + self.L*np.cos(theta1))/(-self.l))
		theta3 = np.arcsin((bvec[0] - x0[0])/(self.l*np.cos(theta2)))
		print(theta1*RAD2DEG,theta2*RAD2DEG, (theta1+theta2)*RAD2DEG)

		#theta3 corresponds to the angle of the ball joint
		if(theta3 > self.swivel_limit):
			return False
		if(np.isnan(theta3) or np.isnan(theta2)):
		    return False
		if(theta1 >= pi/2):
		    return False
		if(theta1 <= -pi/2):
		    return False
		  
		return True


# TEST CODE BELOW
# ds = deltaSolver()
# test = position(0, 0, -0.9)
# print(test.x, test.y, test.z)
# print(ds.ik(test))
# print(ds.vel([1, 1, 1]))

# test2 = position(0.3, 0.5, -1.1)
# print(test2.x, test2.y, test2.z)
# print(ds.ik(test2))

def testPlot():
	kin = deltaSolver()
	kin.plot()
	time.sleep(1)
	kin.updatePlot(position(0, 100, kin.z))
	time.sleep(1)
	kin.updatePlot(position(100, 100, kin.z))
	time.sleep(1)
	kin.updatePlot(position(100, -100, kin.z))
	time.sleep(1)
	kin.updatePlot(position(-100, -100, kin.z))
	time.sleep(1)
	kin.updatePlot(position(-100, 100, kin.z))
	time.sleep(1)
	kin.updatePlot(position(0, 100, kin.z))
	time.sleep(1)
	print(kin.z)
	kin.updatePlot(position(0, 0, kin.z))
	time.sleep(1)

if __name__ == "__main__":
	testPlot()
