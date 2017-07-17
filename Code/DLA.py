import random 
import math



def CirclePoint(cell_list, n=100):
	circle_list=[]
	r_max= cell_list.radius_max
	for i in range(0, n):
		theta= random.uniform(0,1)
		circle_list.append([(r_max+2)*math.cos(theta), (r_max+2)*math.sin(theta)])

	point_coord= random.randint(0,len(circle_list)-1) 
	point= circle_list[point_coord] 
	return point

class cell_list:
	def __init__(self, points):
		self.points=[]
		self.points.append(points) 
		self.radius_max= 1
	def addPoint(self, point):
		self.points.append(point)
		return self.points

	def updateRadiusMax(self):
		max_x= max(self.points, key= lambda x: x[0])[0]
		max_y= max(self.points, key= lambda y: y[1])[1]
		if max_x> max_y:
			self.radius_max= max_x
		else: 
			self.radius_max= max_y
		return self.radius_max

class point:
	def __init__(self, coord):
		self.coord= coord
	def move(self): 
		move_list= [[self.coord[0]-1, self.coord[1]],[self.coord[0]+1, self.coord[1]],[self.coord[0], self.coord[1]-1], [self.coord[0], self.coord[1]+1] ]
		pick_move= random.randint(0, len(move_list)-1)
		self.coord= move_list[pick_move]
		return self.coord

def distance(pt1, pt2):
	return math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)

def isNext(cluster, pot_point):
	for point in cluster:
		if distance(point, pot_point)<=1:
			return True 
	return False 

def getNewPoint(timevalue):
	lst1= cell_list([0,0])
	t=0 
	while t< timevalue:
		temp_point= point(CirclePoint(lst1))
		addable= True 
		while addable and not isNext(lst1.points, temp_point.coord) :
			temp_point.move()
			if temp_point.coord[0]>= lst1.radius_max*3 or temp_point.coord[1]>= lst1.radius_max*3: 
				addable= False 
		if addable: 
			lst1.addPoint(temp_point.coord)
			lst1.updateRadiusMax()
		t+=1  
	return lst1
cells= getNewPoint(800)
f= open("DLA_points.csv", 'w')
f.write("x" + "," + "y" +","+ "frame" +"\n")
count=1
for point in cells.points:
	f.write(str(point[0])+ "," + str(point[1]) + ","+ str(count))
	f.write("\n")
	count+=1
	#ani= animation.FuncAnimation(fig, animate, interval=1000)
f.close()



