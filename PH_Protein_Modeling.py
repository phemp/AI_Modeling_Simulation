import csv
import math
import statistics
import numpy
import random
from random import randint
import matplotlib.pyplot as plt

##############
###Functions
##############


def get_new_pos (prev_pos, step): #Given a position on the lattice and a Manhattan step, calculate the new lattice position
	if step[0] == 0:
		if step[1] == 1:
			##calc new position
			new_pos = [prev_pos[0] , prev_pos[1]+1]
			return new_pos
		
		elif step[1] == 0:
			##calc new position
			new_pos = [prev_pos[0]+1 , prev_pos[1]]
			return new_pos
		
		else:
			##Something has gone wrong
			print ("Manhattan step not correct 1")
				
	elif step[0] == 1:
		if step[1] == 1:
			##calc new position
			new_pos = [prev_pos[0] , prev_pos[1]-1]
			return new_pos
		
		elif step[1] == 0:
			##calc new position
			new_pos = [prev_pos[0]-1 , prev_pos[1]]
			return new_pos
		
		else:
			##Something has gone wrong
			print ("Manhattan step not correct 2")

	else:
		##Something has gone wrong
		print ("Manhattan step not correct 3")

		
def legal (position, d): #Determine if a given lattice position is currently occupied
	#position is not open in the lattice
	if position in d.values():
		return False
	else: 
		return True

		
def folding (iter_start, i, pd, o, iter_end): #folds the sequence
	sto_is = iter_start
	iter = iter_start
	sto_pd = pd.copy()
	p = pd
	op = [[0,0],[0,1],[1,1],[1,0]]
	count = 0
	while iter < iter_end: ##iterate through all the positions that we want folded
		if iter == iter_start: ##This is for if we pass a limited number of options to start with
			options = o[:]
			
		else:
			options = op[:]
			
		check = True
		while (check):
			#Randomly choose Manhattan Step
			if len(options) == 1:
				ri = 0
			elif len(options) == 0:
				print("Something has gone wrong with options")
			else:
				ri = numpy.random.randint(0, len(options)-1)
				
			randstep = options[ri]
			#Fold the amino acid
			fold = get_new_pos(p[iter-1], randstep)
			
			#is the fold legal?
			if legal(fold, p):
				#place the amino acid in the lattice 
				p[iter] = fold
				#Successful fold
				check = False
				
			else:
				#remove move from options
				del options[ri]
				
				#check if we are out of options
				if  len(options) == 0:
					
					#Check if folding is impossible
					if (iter == sto_is):
						return False
						
					#Start all over
					count += 1
					if count == 20:
						return False
					p = sto_pd
					iter = sto_is - 1
					check = False
		
					
		iter += 1
	return p

	
def calc_energy(i, pd, steps):
	energy=0
	for x in range(0, len(i)):
		if i[x] == "H":
			current = pd[x]
			neighbors = []
			
			for y in range (0,len(steps)):
				check_step = steps[y]
				neighbors.append(get_new_pos(current, check_step))
			
			for z in range (0,len(neighbors)):
				#determine if neighbors are occupied
				if neighbors[z] in pd.values():
					#get key that corresponds to the neighbors position
					n_key = list(pd.keys())[list(pd.values()).index(neighbors[z])]
					#check if neighbor is "H"
					if i[n_key] == "H":
						#check if neighbor is topological
						if n_key-1 != x and n_key+1 != x:
							#update energy calc
							energy += -1
	return(energy)


def reconfig (i, pd, steps): #Outputs new configuration 
	pop_step = steps[:]
	#eliminate the step used in the previous config, so that we are guaranteed to propose new change
	for s in range(0,len(steps)):
		if pd[len(i)-1] == get_new_pos(pd[(len(i)-1)], pop_step[s]) :
			del pop_step[s]
	
	
	npd = folding(len(i)-1, i, pd, pop_step, len(pd))
	return npd	
	
def point_mut (i, pd): ##create a point mutation and outputs the resulting configuration
	count_ri = 0
	step = [[0,0],[0,1],[1,1],[1,0]]
	##Choose random index to change##
	randchange = randint(2,len(i)-2)
	count_ri += 1
	new_c = []

	##Get old configuration up to where the change will be made##
	for y in range(0,len(i)):
		if y < randchange:
			new_c.append(i[y])

	##Get new configuration##
	new_posdict = reconfig(new_c, pd, step)
	
	#If the fold was impossible
	while new_posdict == False:
		if count_ri == (len(i)-2):
			return False
		#Re-do
		randchange = randint(2,len(i)-2)
		count_ri +=1
		new_c = []
		
		for y in range(0,len(i)):
			if y < randchange:
				new_c.append(i[y])
				
		new_posdict = reconfig(new_c, pd, step)
		
	return new_posdict

def simulate(iter_end, i, pd, c, T): #Runs a number of simulations and outputs the energies in a list
	sample_energy = []
	old_pd = pd
	old_config_energy = c
	for x in range(0,iter_end):
		new_pd = point_mut(i, old_pd)
		step = [[0,0],[0,1],[1,1],[1,0]]
		
		
		##Calc Energy##
		new_config_energy = calc_energy(i, new_pd, step)
		
		##Calc Metropolis Alg Values##
		n_nce = -1 * new_config_energy  #Neg new energy
		n_oce = -1 * old_config_energy  #Neg old energy 
		
		met_value = (math.exp(n_nce/T) / math.exp(n_oce/T))
		
		##Metropolis accept/reject##
		if new_config_energy < old_config_energy:
			#ACCEPT
			sample_energy.append(new_config_energy)
			
			
			##Set old_pd to new_pd for the next iteration
			old_pd = new_pd
			old_config_energy = new_config_energy
		
		elif random.uniform(0,1) < met_value:
			#ACCEPT
			sample_energy.append(new_config_energy)
			
		
			##Set old_pd to new_pd for the next iteration
			old_pd = new_pd
			old_config_energy = new_config_energy
			
		else:
			##REJECT
			sample_energy.append(old_config_energy)
	
	return sample_energy	
	

########
##MAIN##
########
	
	
#############################
###Read in the file to a list
#############################

file = open("hw1.data.txt", "rU") 
input = []
for line in file:
	for c in line:
		input.append(c)

#testable sim input
#input = ["H","H","P","P","H","H","P","P"]


###################################################
###Part 1 - Build the lattice model and calc energy
###################################################

start_dict = {0 : [0,0]} #Holds lattice positions
man_steps = [[0,0],[0,1],[1,1],[1,0]]

##Fold the next amino acid in the input sequence##

position_dict = folding(1, input, start_dict, man_steps, len(input))
	
##See output for debugging
#for i in position_dict:
	#print (i, position_dict[i])
	
##Calc Energy for model##
config_energy = calc_energy(input, position_dict, man_steps)
print(config_energy)



#############################################						
###Part 2 - 1000 simulations of the model 
#############################################


s = simulate(1000, input, position_dict, config_energy, 100)


plt.hist(s, align='mid', bins=range(min(s), max(s) + 2, 1))
plt.title("Histogram of Energy: 1000 Configurations")
plt.xlabel("Energy")
plt.ylabel("Count")
plt.show()


#########################
###Part 3 - Optimization 
#########################

#Initialize Population
population = [] #A list of configurations (dictionaries)
best_energy = 0 #Inititalize best energy variable
best_config = {0 : [0,0]}

for x in range(0,100): #Population the population and get the value of the best energy
	pop_f = folding(1, input, start_dict, man_steps, len(input)) 
	if pop_f == False:
		x = x-1
		break
	
	else:
		population.append(pop_f)
	
	ce = calc_energy(input, population[x], man_steps)
	if ce <= best_energy:
		best_energy = ce
		best_config = population[x]
	
		
#Optimization		
t = 70
current_pop = population
next_gen = []
check = True
while (check):
	n=0
	y = randint(0,len(current_pop)-1)
	#Make a pointwise mutation of the next two configurations
	pm_config_m = current_pop[y] 
	
	pm_config_f = point_mut(input, current_pop[y])
	
	if(pm_config_f == False):
		y = randint(0,len(current_pop)-1)
		pm_config_m = current_pop[y] 
		pm_config_f = point_mut(input, current_pop[y])
	
	
	while (n < len(current_pop)): #The second generation should make an equal or less number of individuals in the initial population
		Produce child via crossover
		cross_loc = randint(1, len(input)-2)
		
		
		child = {0 : [0,0]}
		for z in range(1,len(input)):
			if z <= cross_loc:
				child[z] = pm_config_m[z]
			
			else:
				if legal(pm_config_f[z], child):
					child[z] = pm_config_f[z]
				
				else: #Crossover is not legal
					break
		
		if len(child) == len(pm_config_m): ##We sucessfully made a legal child
			child_e = calc_energy(input, child, man_steps)
		
			m_e = calc_energy(input, pm_config_m, man_steps)
			f_e = calc_energy(input, pm_config_f, man_steps)

			p_av = statistics.mean([m_e, f_e])
			
			if (child_e <= p_av):
				next_gen.append(child)
				
				
			else:
				met_rand = random.uniform(0,1)
				if (met_rand < (-(p_av + child_e))/t):
					next_gen.append(child)
		
			
			if (child_e <= best_energy):
				best_energy = child_e
				best_config = child
		
		n += 1
	t = t/math.log(n)
	
	if len(current_pop) < 2 or t == 0:
		check = False
		break
	
	
	#Update the generations and prepare to make another
	current_pop = next_gen
	next_gen = []
	
	if len(current_pop) < 2 or t == 0:
		check = False

##Write best config to a file
with open("dict.csv", "w") as csv_file:
	writer = csv.writer(csv_file)
	for key, value in best_config.items():
		writer.writerow([key, value])


print ("BEST ENERGY")
print (best_energy)


file.close()