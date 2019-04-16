import numpy as np
from numpy import arctan, arctan2, sin, cos, tan, sqrt

"""
Solve the inverse kinematics to get desired theta1, theta2, and theta3

INPUT:	goal: [xe, ye, ze], current_thetas: [theta1, theta2, theta3]
OUTPUT:  desired_thetas: [theta1, theta2, theta3]

Length of links = [a1, a2, a3]

[B1, B2, B3] = Mounting points of servos on base plate
[P1, P2, P3] = Tips of link arms on end effector plate
E = center of end effector plate

3-dimensional coordinates of these vectors:
OB = [OB1, OB2, OB3]
BA = [B1A1, B2A2, B3A3]
EP = [EP1, EP2, EP3]
AP = [A1P1, A2P2, A3P3]
OE = [xe, ye, ze]

L = length of upper arm
l = length of lower arm (the parallelogram)
s = edge length of base plate triangle
d = edge length of end effector triangle
E = [xe, ye, ze] = center of end effector plate
O = origin = center of base plate

3 kinematic conditions based on:
OE = OB + BA + AP + PE
AP = l

in the form of:
a * cos(theta) + b * sin(theta) + c = 0
a = [a1, a2, a3]
b = [b1, b2, b3]
c = [c1, c2, c3]
theta = [theta1, theta2, theta3]

Use tangent half-angle formula to rewrite kinematic constraint eqns as 2nd-order algebraic eqns

Obtain conditions under which the IK problem has a soln.
"""

s = 20	# Base plate triangle edge length
d = 5	# End effector triangle edge length
L = 30	# Upper arm length
l = 25	# Lower arm length (parallelogram)

def in_joint_range(q):
    for i, qi in enumerate(q):
        if qi < joint_limits[i][0] or qi > joint_limits[i][1]:
            return False
    return True

def select_best_candidate(candidates, q0, weight = [1, 1]):
	# we prefer minimum travel
    min_v = None
    best_q = None
    for q in candidates:
        v = np.sum(((np.array(q) - np.array(q0)) * np.array(weight)) ** 2)
        if (min_v == None or min_v > v) and in_joint_range(q):
            min_v = v
            best_q = q
    return best_q

def solve_inverse_kinematics(goal, current_thetas):
	OE = goal
	theta1 = current_thetas[0]
	theta2 = current_thetas[1]
	theta3 = current_thetas[2]

	OB1 = [0, -s/(2*sqrt(3)), 0]
	B1A1 = [0, -L * cos(theta1), -L * sin(theta1)]
	EP1 = [0, -d/sqrt(3), 0]

	OB2 = [s/4, s/(4*sqrt(3)), 0]
	B2A2 = [sqrt(3)/2 * L * cos(theta2), 1/2 * L * cos(theta2), -L * sin(theta2)]
	EP2 = [d/2, d/(2*sqrt(3)), 0]

	OB3 = [-s/4, s/(4*sqrt(3)), 0]
	B3A3 = [-sqrt(3)/2 * L * cos(theta3), 1/2 * L * cos(theta3), -L * sin(theta3)]
	EP3 = [-d/2, d/(2*sqrt(3)), 0]

	OB = [OB1, OB2, OB3]
	BA = [B1A1, B2A2, B3A3]
	EP = [EP1, EP2, EP3]

	ax, ay, bx, by, px, py, X, Y, a, b, c = ([] for i in range(11))
	for i in range(0, 3):
		ax.append(BA[i][0])
		ay.append(BA[i][1])
		bx.append(OB[i][0])
		by.append(OB[i][1])
		px.append(EP[i][0])
		py.append(EP[i][1])
		X.append(bx[i] + px[i] - OE[0])
		Y.append(by[i] + py[i] - OE[1])
		a.append(2 * (ax[i] * X[i] + ay[i] * Y[i]))
		b.append(2 * L * OE[2])
		c.append(X[i]**2 + Y[i]**2 + OE[2]**2 - l**2 + L**2)

		if (a[i]**2 + b[i]**2 - c[i]**2 >= 0):	# For all i, t must be real
			candidates = [(-b[i] + sqrt(a[i]**2 + b[i]**2 - c[i]**2)) / (c[i] - a[i]),
						  (-b[i] - sqrt(a[i]**2 + b[i]**2 - c[i]**2)) / (c[i] - a[i])]
			t[i] = select_best_candidate(candidates, current_thetas)
		else:
			t[i] = None

		desired_thetas[i] = 2 * arctan(t[i])

	return desired_thetas

goal = [10, 10, 10]
current_thetas = [1/sqrt(2), 1/sqrt(2), 1/sqrt(2)]
solve_inverse_kinematics(goal, current_thetas)









