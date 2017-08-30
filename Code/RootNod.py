
#create root nodule colonization areas
class Place(object): 
	def __init__(self, x, y):
		self.x_coord=x 
		self.y_coord=y 
		self.coordinate= [x,y]
		self.tenant= None 
		self.hasBacteria=False

	def add_bacteria(self, bacteria):
		if self.tenant is None: 
			self.tenant= bacteria
			bacteria.place=self
			self.hasBacteria=True



def CirclePoint(cell_list, n=100):
	circle_list=[]
	r_max= cell_list.radius_max
	for i in range(0, n):
		theta= random.uniform(0,1)*2*math.pi
		circle_list.append([(r_max+2)*math.cos(theta), (r_max+2)*math.sin(theta)])

	point_coord= random.randint(0,len(circle_list)-1) 
	coordinate= circle_list[point_coord] 
	return coordinate

class particle:
	def __init__(self, coord):
		self.time_alive=0 
		self.coord= coord
		
	def move(self): 
		move_list= [[self.coord[0]-1, self.coord[1]],[self.coord[0]+1, self.coord[1]],[self.coord[0], self.coord[1]-1], [self.coord[0], self.coord[1]+1] ]
		pick_move= random.randint(0, len(move_list)-1)
		self.coord= move_list[pick_move]
		return self.coord

def distance(pt1, pt2):
	return math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)

#checks if aggregateis next to particle
def isNext(cluster, pot_point):
	for point in cluster:
		if distance(point.coord, pot_point)<=1:
			return True 
	return False 

class cell_list:
	def __init__(self, points):
		self.points=[]
		self.points.append(points) 
		self.radius_max= 1
	def addCell(self, cell):
		self.points.append(cell.coord)
		return self.points
	def updateRadiusMax(self):
		max_x= max(self.points, key= lambda x: x.coord[0]).coord[0]
		max_y= max(self.points, key= lambda y: y.coord[1]).coord[1]

		if max_x> max_y:
			self.radius_max= max_x
		else: 
			self.radius_max= max_y
		return self.radius_max

class Bacteria(object):
	def __init__(self, growth_type, place=None, color):
		self.growth_type= growth_type
		self.place=None 
		self.color=color
	def getNewPoint(timevalue):


class FilamentousBacteria(Bacteria): 
	def getNewPoint(self, timevalue):
		cell_lst= cell_lst(self.coord) 
		t=0 
		while t<timevalue:
			temp_particle= particle(CirclePoint(cell_lst))
			addable= True 
			while addable and not isNext(cell_lst.points, temp_particle.coord): 
				temp_point.move()
				if temp_particle.coord[0]>= lst1.radius_max*3 or temp_particle.coord[1]>= lst1.radius_max*3: 
					addable= False
			if addable: 
				new_fil_cell= FilamentousBacteria("filamentous", red)
				cell_lst.addCell(temp_particle)
				place_list[temp_particle.coord].add_bacteria(new_fil_cell)



place_list=[]
for x in range(0, 10):
	for y in range(0, 10):
		place_list.append(Place(x,y))

new_bac= Bacteria("filamentous", red)
print(new_bac.place)
place_list[1].add_bacteria(new_bac)