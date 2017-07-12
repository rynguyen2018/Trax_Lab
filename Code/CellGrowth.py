from random import randint
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib import style 

style.use('fivethirtyeight')



point= [0,0]
black_list= []
black_list.append(point)

def animate(i):
	graph_data=open('Cell_growth_points.csv', 'r').read()
	lines = graph_data.split('\n')
	xs=[]
	ys=[]
	for line in lines: 
		if len(line)>1: 
			x,y= line.split(',')
			xs.append(x)
			ys.append(y)
	ax1.clear()
	ax1.plot(xs, ys)

def getAdjacentPoints(starting_point):
	point1= [starting_point[0]-1, starting_point[1]] 
	point2= [starting_point[0]+1, starting_point[1]]  
	point3= [starting_point[0], starting_point[1]-1]  
	point4= [starting_point[0], starting_point[1]+1]

	return [point1, point2, point3, point4]

def GardenOfEden(black_list, t): 
	print("propagating points")
	while t>0: 
		starting_point= black_list[randint(0, len(black_list)-1)]  
		point_list= getAdjacentPoints(starting_point)
		potential_black_point_list=[]
		for point in point_list: 
			if point not in black_list: 
				potential_black_point_list.append(point)
		if len(potential_black_point_list)!=0: 
			pick_point= potential_black_point_list[randint(0, len(potential_black_point_list)-1)]
			black_list.append(pick_point)
			t-=1 
	return black_list

black_list= GardenOfEden(black_list, 5000)
f= open("Cell_growth_points.csv", 'w')
f.write("x" + "," + "y" +","+ "frame" +"\n")
count=1
for point in black_list:
	f.write(str(point[0])+ "," + str(point[1]) + ","+ str(count))
	f.write("\n")
	count+=1
	#ani= animation.FuncAnimation(fig, animate, interval=1000)
f.close()

# fig= plt.figure()
# ax1= fig.add_subplot(1,1,1)
# ani= animation.FuncAnimation(fig, animate, interval=1000)
# plt.show()


#GardenOfEden(point, black_list, 100)