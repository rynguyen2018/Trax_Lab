import random 
import math
import sys

def CirclePoint(cell_list, n=100):
	circle_list=[]
	r_max= cell_list.radius_max
	for i in range(0, n):
		theta= random.uniform(0,1)*2*math.pi
		circle_list.append([(r_max+2)*math.cos(theta), (r_max+2)*math.sin(theta)])

	point_coord= random.randint(0,len(circle_list)-1) 
	coordinate= circle_list[point_coord] 
	return coordinate

class cell_list:
	def __init__(self, points):
		self.points=[]
		self.points.append(points) 
		self.radius_max= 1
	def addPoint(self, point):
		self.points.append(point)
		return self.points

	def updateRadiusMax(self):
		max_x= max(self.points, key= lambda x: x.coord[0]).coord[0]
		max_y= max(self.points, key= lambda y: y.coord[1]).coord[1]

		if max_x> max_y:
			self.radius_max= max_x
		else: 
			self.radius_max= max_y
		return self.radius_max

class pointer:
	def __init__(self, coord):
		self.time_alive=0 
		self.coord= coord
		
	def move(self): 
		move_list= [[self.coord[0]-1, self.coord[1]],[self.coord[0]+1, self.coord[1]],[self.coord[0], self.coord[1]-1], [self.coord[0], self.coord[1]+1] ]
		pick_move= random.randint(0, len(move_list)-1)
		self.coord= move_list[pick_move]
		return self.coord
	def updateTime(self):
		self.time_alive+=1
		return self.time_alive

def distance(pt1, pt2):
	return math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)

def isNext(cluster, pot_point):
	for point in cluster:
		if distance(point.coord, pot_point)<=1:
			return True 
	return False 

def getAdjacentPoints(starting_point):
	point1= [starting_point[0]-1, starting_point[1]] 
	point2= [starting_point[0]+1, starting_point[1]]  
	point3= [starting_point[0], starting_point[1]-1]  
	point4= [starting_point[0], starting_point[1]+1]

	return [point1, point2, point3, point4]

def Budding(point, cell_list):
	point_list= getAdjacentPoints(point.coord)
	black_list=[]
	potential_black_point_list=[]
	for point in cell_list.points: 
		black_list.append(point.coord)
	for point in point_list:
		if point not in black_list:
			potential_black_point_list.append(pointer(point))
	if len(potential_black_point_list)!=0: 
		#print("new point on the block")
		pick_point= potential_black_point_list[random.randint(0, len(potential_black_point_list)-1)]
	else:
		pick_point= pointer([0,0])
	return pick_point		


def getNewPoint(timevalue):
	lst1= cell_list(pointer([0,0]))
	t=0 
	while t< timevalue:
		temp_point= pointer(CirclePoint(lst1))
		addable= True 
		while addable and not isNext(lst1.points, temp_point.coord) :
			temp_point.move()
			if temp_point.coord[0]>= lst1.radius_max*3 or temp_point.coord[1]>= lst1.radius_max*3: 
				addable= False 
		if addable: 
			lst1.addPoint(pointer(temp_point.coord))
			lst1.updateRadiusMax()
			t+=1 
			print(t)
		for point in lst1.points: 
			point.updateTime()
			if point.time_alive>=10: 
				lst1.addPoint(Budding(point, lst1))
				point.time_alive=0

	return lst1


print("now getting cells")
cells= getNewPoint(30)





f= open("DLA_points.csv", 'w')
f.write("x" + "," + "y" +","+ "frame" +"\n")
count=1
for point in cells.points:
	f.write(str(point.coord[0])+ "," + str(point.coord[1]) + ","+ str(count))
	f.write("\n")
	count+=1
	#ani= animation.FuncAnimation(fig, animate, interval=1000)
f.close()



