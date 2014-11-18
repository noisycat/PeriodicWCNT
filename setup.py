#!/usr/bin/env python
from sympy import Matrix, sqrt, sin, cos, pi, oo, atan, Rational
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
def rotate(coords, theta,_about = (0,0), rounding=False):
	cosa = cos(theta);
	sina = sin(theta)
	xf =  cosa * (coords[0]-_about[0]) - sina * (coords[1]-_about[1]) + _about[0];
	yf =  sina * (coords[0]-_about[0]) + cosa * (coords[1]-_about[1]) + _about[1];
	if rounding:
		xf = round(xf)
		yf = round(yf)
	return Matrix([xf, yf, coords[2]])

def dist(pa,pb):
	return sum([(pa[i]-pb[i])**2 for i in xrange(0,3)])

def NNCalculation(points):
	lsegs = list()
	mindist = min([dist(points[0],points[j]) for j in range(1,len(points))])
	for i in xrange(len(points)):
		l = [(j,abs(dist(points[i],points[j])-mindist) <= 2e-10) for j in xrange(i+1,len(points))]
		for entry in l:
			if entry[1]: 
				lsegs.append(LineSeg(points[i],points[entry[0]]))
	return lsegs

def intersect(l1_coords,l2_coords):
	y_2_1 = l1_coords[1][1]-l1_coords[0][1]
	x_2_1 = l1_coords[1][0]-l1_coords[0][0]
	x_4_3 = l2_coords[1][0]-l2_coords[0][0]
	y_4_3 = l2_coords[1][1]-l2_coords[0][1]
	x_3_1 = l2_coords[0][0]-l1_coords[0][0]
	y_3_1 = l2_coords[0][1]-l1_coords[0][1]
	den_det = (x_4_3 * y_2_1 - x_2_1 * y_4_3)
	if den_det==0: 
		# check to see if same line *above test fails because of nullspace i think*
		try: 
			s = x_3_1/x_2_1
			if y_2_1*s == y_3_1 and 1 >= s >= 0: return (s,numpy.inf,True)
		except: 
			s = y_3_1/y_2_1
			if x_2_1*s == x_3_1 and 1 >= s >= 0: return (s,numpy.inf,True)
		try:
			s = (l2_coords[1][0]-l1_coords[0][0])/x_2_1
			if y_2_1*s == (l2_coords[1][1]-l1_coords[0][1]) and 1 >= s >= 0: return (s,numpy.inf,True)
		except: 
			s = (l2_coords[1][1] - l1_coords[0][1])/y_2_1
			if x_2_1*s == (l2_coords[1][0]-l1_coords[0][0]) and 1 >= s >= 0: return (s,numpy.inf,True)

		# if not same line, then parallel dont touch awayyyyy
		return (-oo,-oo,False)
	s = (x_4_3 * y_3_1 - x_3_1 * y_4_3) /  den_det
	t = (x_2_1 * y_3_1 - x_3_1 * y_2_1) /  den_det
	flag = True if (s <= 1 and s >= 0 and t <= 1 and t >=0) else False
	if flag: print (float(s),float(t),flag)
	return (float(s),float(t),flag)

class Atom:
	def __init__(self, point, name):
		self.coords = point
		self.atom   = name

	def __getitem__(self,i):
		return self.coords[i]

class Ring:
	def __init__(self, center, rotset = 0):
		v = Matrix([1,0,0])
		self.atoms = [center + a/2.0 * rotate(v,0+rotset),
			center + a/2.0 * rotate(v,pi/3+rotset),
			center + a/2.0 * rotate(v,2*pi/3+rotset),
			center + a/2.0 * rotate(v,pi+rotset),
			center + a/2.0 * rotate(v,4*pi/3+rotset),
			center + a/2.0 * rotate(v,5*pi/3+rotset)]

	def Translate(self,x,y):
		m = Matrix([x,y,0])
		for i in xrange(len(self.atoms)):
			self.atoms[i] += m

class LineSeg:
	def __init__(self,point1,point2):
		self.coords = (point1,point2)
	
	def __getitem__(self,i):
		return self.coords[i]

	def plot(self,arrow=False):
		if arrow:
			plt.arrow(float(self.coords[0][0]),float(self.coords[0][1]),float(self.coords[1][0]-self.coords[0][0]),float(self.coords[1][1]-self.coords[0][1]),'-g')
		else:
			plt.plot([self.coords[0][0],self.coords[1][0]],[self.coords[0][1],self.coords[1][1]],'-g')

	def plot3(self):
		plt.plot([self.coords[0][0],self.coords[1][0]],[self.coords[0][1],self.coords[1][1]],[self.coords[0][2],self.coords[1][2]],'-g')

	def Intersect(self,LineSeg_):
		return intersect(self.coords, LineSeg_.coords)

	def Rotate(self,theta,about = (0,0,0)):
		self.coords = (rotate(self.coords[0], theta, about), rotate(self.coords[1],theta, about))
	
class LineLoop:
	def __init__(self,*points):
		self.coords = points

	def plot(self):
		x = [l[0] for l in self.coords]; x.append(x[0])
		y = [l[1] for l in self.coords]; y.append(y[0])
		plt.plot(x,y,'--g')

	def Contains(self,point):
		intersections = 0
		# construct ray from point
		ray = LineSeg(point,Matrix([10000,1,0]))
		# test intersection of ray with all sides.
		for i in range(len(self.coords)-1):
			s,t,flag = ray.Intersect(LineSeg(self.coords[i],self.coords[i+1]))
			if s < 2e-16 and flag: return True
			intersections += 1 if flag else 0
		# closing loop
		s,t,flag = ray.Intersect(LineSeg(self.coords[-1],self.coords[0]))
		if s < 2e-16 and flag: return True
		intersections += 1 if flag else 0
		# count. if odd, inside, if even, outside
		return (intersections % 2 == 1)

def transform(r,theta,l): 
	x = r * cos(theta)
	y = r * sin(theta)
	z = l
	return [x,y,z]

a = Rational(2.46) #angstrom
Av = [ LineSeg(Matrix((0,0,0)),Matrix([sqrt(3)*a/2, a/2, 0])), LineSeg(Matrix((0,0,0)),Matrix([sqrt(3)*a/2, -a/2, 0])) ]
b = [ Matrix([0,0,0]),Matrix([a/(2*sqrt(3)),a/2,0])]
theta = Rational(0)

A = [Av[0].coords[1], Av[1].coords[1]]
b = [ Matrix([0,0,0]),Matrix(rotate([a/(2*sqrt(3)),a/2,0],theta))]
def ConstructSWNT(n = 6, m = 0, show=False, Repeated=Rational(1), DebugPoints=False, offset=Matrix([0,0,1]),rotY=Rational(0)):
	R = list()
	lsegs = list()
	Ri = Rational(n) * A[0] + Rational(m) * A[1]
	Ch = LineSeg(Matrix((Rational(0),Rational(0),Rational(0))),Ri)
	def _gcd(a,b): 
		'''Calculate greatest common divisor of a,b'''
		while b: 
			a, b = b, a % b
		return a
	d = _gcd(n, m)
	d_r = 3*d if (((n-m) % (3*d)) == 0) else d
	# anything at all that seems special comes from flex.phys.tohoku.ac.jp/~rsaito/rsj/f578.pdf - "Physics of Carbon Nanotubes"
	t = [ (2*m+n) / d_r , -(2*n+m)/d_r]
	print d,d_r, t,Matrix(t)+Matrix([n,m])
	T = LineSeg(Matrix([0,0,0]),t[0] * A[0] + t[1] * A[1])
	print T.coords[1]

	tolm = Matrix([Rational(0),Rational(0),Rational(0)])
	Domain = LineLoop(Ch.coords[0]-tolm, Ch.coords[1]+tolm, T.coords[1]+Ch.coords[1]+tolm, T.coords[1]+Ch.coords[0]+tolm)

	theta = atan(Ch.coords[1][1]/Ch.coords[1][0])
	plt.figure(0).add_subplot(211)
	plt.title("Slice of coordinates")
	nlimits = (min([t[0],0]),max([t[0]+n]))
	mlimits = (min([t[1],0]),max([t[0]+m]))
	for i in xrange(nlimits[0],nlimits[1]+1):
		for j in xrange(mlimits[0],mlimits[1]+1):
			for bi in b: 
				R.append(bi + i * A[0] + j * A[1])
				if Domain.Contains(R[-1]) and DebugPoints: lsegs.extend([LineSeg(bi,bi+i*A[0]),LineSeg(bi+i*A[0],bi+i*A[0]+j*A[1])])

	plt.plot([x[0] for x in b], [x[1] for x in b],'rx')
	plt.plot([x[0] for x in R], [x[1] for x in R],'ob')
	Rselected = [x for x in R if Domain.Contains(x)]
	plt.plot([x[0] for x in Rselected], [x[1] for x in Rselected],'og')
	Ch.plot(arrow=True)
	T.plot(arrow=True)
	Av[0].plot(arrow=True)
	Av[1].plot(arrow=True)
	Domain.plot()
	plt.axis('equal')

	#lsegs = NNCalculation(Rselected)
	for lseg in lsegs: lseg.plot()
#rotate to Ch to xaxis
	theta = atan(Ch.coords[1][1]/Ch.coords[1][0])
	T.Rotate(-theta)
	Ch.Rotate(-theta)
	for i in xrange(len(Rselected)): 
		Rselected[i] = rotate(Rselected[i],-theta)
	Rselected = [x+T.coords[1]*i for x in Rselected for i in xrange(Repeated)]
	plt.figure(0).add_subplot(212)
	plt.title('Selected After Rotation transform')
	for lseg in lsegs:
		lseg.Rotate(-theta)
		lseg.plot()
	plt.plot([x[0] for x in Rselected], [x[1] for x in Rselected],'og')
	plt.axis('equal')
	fig = plt.figure(2)
	ax = fig.add_subplot(211,projection='3d')
	ax.plot( [x[0] for x in Rselected], [x[1] for x in Rselected], [x[2] for x in Rselected], 'o')
	x0 = Rational(0)
	xf = Ch.coords[1][0]
	def thetaf(x): return pi * (Rational(2)*(x-x0)/(xf-x0)-Rational(1))
	def rf(x): return (xf-x0)/(Rational(2)*pi)
	#def xx(x): return rf(x) * (cos(thetaf(x)-pi/2)) + (xf+x0)/2
	def xx(x): return rf(x) * (cos(thetaf(x)-pi/2)) + offset[0]
	def yy(x): return rf(x) * (sin(thetaf(x)-pi/2)  +  offset[2])
	def transformR(x): return Matrix([xx(x[0]), x[1], yy(x[0])])

	ax.plot( [xx(x[0]) for x in Rselected], [x[1] for x in Rselected], [yy(x[0]) for x in Rselected], 'o')
	RTransformed = list()
	for x in Rselected:
		tmprot = rotate((xx(x[0]),yy(x[0]),0),rotY)
		RTransformed.append(Matrix([tmprot[0], x[1], tmprot[1]]))

	for lseg in lsegs: 
		lseg.coords = [[xx(lseg.coords[i][0]) , lseg.coords[i][1], yy(lseg.coords[i][0])] for i in range(2)]
		lseg.plot3()

	ax.axis('equal')
	plt.title('Transformation into SWNT')
	ax = fig.add_subplot(212,projection='3d')
	ax.plot( [x[0] for x in Rselected], [x[1] for x in Rselected], [x[2] for x in Rselected], 'o')
	x0 = 0
	xf = Ch.coords[1][0]
	ax.plot( [xx(x[0]) for x in Rselected], [x[1] for x in Rselected], [yy(x[0]) for x in Rselected], 'o')
	for lseg in lsegs: 
		lseg.plot3()
	ax.axis('equal')
	plt.title('Transformation into SWNT')
	from os import mkdir, path
	foldername = 'n,m-%d,%d' % (n,m)
	try: mkdir(foldername)
	except: pass
	plt.figure(0).savefig(path.join(foldername,'graphene-%d,%d.png' % (n,m)))
	plt.figure(2).savefig(path.join(foldername,'tube-%d,%d.png' % (n,m)))
	if show: plt.show()
	

	#RTransformed = [transformR(x) for x in Rselected]
# find duplicate points after transformation and remove them
	RDuplicates = list()
	for i in range(len(RTransformed)-1,0,-1):
		if RTransformed.count(RTransformed[i]) > 1: RDuplicates.append(RTransformed.pop(i))
	num_tube_duplicates = len(RDuplicates)
# find points that are duplicates in the unit cell, and remove them for the file print,
	Rtranslated = [x+Repeated*T[1] for x in RTransformed if x+Repeated*T[1] in RTransformed]
	num_tube_duplicates = num_tube_duplicates + len(Rtranslated)
	print 'points removed #:',num_tube_duplicates
	for xT in Rtranslated: RTransformed.remove(xT)
	plt.figure(0).clf()
	plt.figure(2).clf()
	return RTransformed, T, xx, yy

# rearrange these so that the lattice is in the x direction for repeating units
def tempOut(Rselected, n, m, T, xx=lambda i: i, yy=lambda i: i, foldername='',order=(0,1,2)):
	from os import path,mkdir
	foldername= 'n,m-%d,%d' % (n,m) if foldername == '' else foldername
	try: mkdir(foldername)
	except: pass
	parsecin = open(path.join(foldername,'parsec.in'),'w')
	parsecin.write('''
Boundary_Conditions wire

Boundary_Sphere_Radius 10 ang

begin Cell_shape
%24.16lf %24.16lf %24.16lf
end cell_shape

lattice_vector_scale 1 ang

Grid_Spacing:  0.55

States_Num:  %d

eigensolver chebdav

coordinate_unit Cartesian_Ang

Correlation_Type ca

Ignore_Symmetry true

Max_Iter 100

Atom_Types_Num:  1

#------------ new atom type ----------------
Atom_Type: C
Local_Component: s

begin Atom_Coord\n''' % (float(T[1][order[0]]),float(T[1][order[1]]),float(T[1][order[2]]),len(Rselected)*6))
	for x in Rselected: 
		parsecin.write('%24.16lf %24.16lf %24.16lf\n' % (float(x[order[0]]),float(xx(x[order[1]])),float(yy(x[order[2]]))))
	parsecin.write('end Atom_Coord')
	parsecin.close()
	xcrysden = open(path.join(foldername,'out.xyz'),'w')
	xcrysden.write('\t%d\n\n' % (len(Rselected)))
	for x in Rselected: 
		xcrysden.write('C %24.16lf %24.16lf %24.16lf\n' % (float(x[order[0]]),float(xx(x[order[1]])),float(yy(x[order[2]]))))
	xcrysden.close()

def RotateSWNT(coords,theta):
	RTransformed = list()
	for x in coords:
		tmprot = rotate((x[0],x[2],x[1]),theta,_about=(0,0))
		RTransformed.append(Matrix([tmprot[0], x[1], tmprot[1]]))
	return RTransformed

def MWNT_Rots(nmcollection, rotationdict, Repeated=Rational(1)):
	Tube = list()
	i = 0
	for n,m in nmcollection:
		iTube,T,_,_ = ConstructSWNT(n,m,Repeated=Repeated,show=False,offset=Matrix((0,0,0)))
		Tube.append(RotateSWNT(iTube, rotationdict.get((n,m),Rational(0))))
		i = i+1
	Tubes = list(itertools.chain.from_iterable([R for R in Tube]))
	return Tubes,Tube,T

class MWNT:
	def __init__(self,nmcollection,rotation=dict(),Repeated=Rational(1)):
# nm collection is a list of (n,m) tuples
# rotationdict is a keylist of (n,m) : theta pairs
		self.T = dict()
		self.Tube = dict()
		self.Theta = dict()
		for n,m in nmcollection:
			iTube,T,_,_ = ConstructSWNT(n,m,Repeated=Repeated,show=False,offset=Matrix((0,0,0)))
			self.Tube[(n,m)] = iTube
			self.T[(n,m)] = T
			self.Theta[(n,m)] = rotation.get((n,m),0)

	def Flatten(self):
		Tubes = list(itertools.chain.from_iterable(self.Tube.values()))
		return Tubes

	def Rotate(self, rotation):
		for n,m in rotation:
			self.Tube[(n,m)] = RotateSWNT(self.Tube[(n,m)], rotation.get((n,m), Rational(0)))

	def Export(self):
		# temp way to do this
		maxn, maxm = self.Tube.keys()[-1]
		outstr = 'nm-'+','.join(('%dx%d' % (n,m) for (n,m) in self.Tube.keys())) + '-%f' % (0)
		tempOut(self.Flatten(),maxn,maxm,self.T[(maxn,maxm)], foldername=outstr, order = (1,0,2))



#for k in [3,5,8]:
#	Construct = MWNT([(k,0)],{(k,0):0})
#	tempOut(Construct.Flatten(), k, 0, Construct.T[(k,0)], foldername='nm-%dx%d-%f' % (k,0,0), order=(1,0,2))
#	rot = 2*pi/(k*10)
#	for i in xrange(1,6):
#		Construct.Rotate({(k,0):rot})
#		tempOut(Construct.Flatten(), k, 0, Construct.T[(k,0)], foldername='nm-%dx%d-%f' % (k,0,rot*i), order=(1,0,2))

#for i in xrange(0,5):
#	rot = 2*pi/20*i
#	Construct,T = MWNT([(8,2),(16,2)],{(8,2):rot})
#	tempOut(Construct,16,2,T,foldername='nm-8,16x2-%f' % (rot))


#for i in range(3,9):
	#for j in range(0,9): ConstructSWNT(i,j)

# testing d, d_r
#for i,j in [(5,5),(9,0),(6,5),(7,4),(8,3),(10,0)]:
	#ConstructSWNT(i,j,show=True,Repeated=1,DebugPoints=True)
#apply transformations to wrap around
# (5,0) x (10,0)
Construct = MWNT([(5,0),(10,0)],Repeated=Rational(2))
Construct.Export()

Construct = MWNT([(8,0),(17,0)])
Construct.Export()

Construct = MWNT([(8,0),(16,0)])
Construct.Export()
# (100 atom system)
# (1000 atom system)
# (5,0)
#rot = 0
#Construct,T = MWNT([(5,0),(10,0)],{(5,0):0})
#tempOut(Construct,10,0,T,foldername='nm-5,10x0-%f' % (rot))
## (8,0)
#Construct,T = MWNT([(5,0)],{(5,0):0})
#tempOut(Construct,5,0,T,foldername='nm-5x0-%f' % (rot))
## (10,0)
#Construct,T = MWNT([(10,0)],{(5,0):0})
#tempOut(Construct,10,0,T,foldername='nm-10x0-%f' % (rot),order=(1,0,2))
#
