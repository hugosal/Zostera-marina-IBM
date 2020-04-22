#o-*- coding: utf-8 -*-
#Hugo Salinas 2019

#function to plot a simulation
from zostera_plotter import plot_meadow

import random as rn
import os
import sys
import pickle
import pandas
import math
import csv
import numpy as np
from datetime import datetime
import time
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import gamma
from scipy.stats import poisson
from scipy.stats import bernoulli
from tqdm import tqdm # for progress bar for loop
# the zostera class represents a individual of z. marina. consisting of at least one branch
class zostera:
	def __init__(self, branches):  # branches is a list of branch class instances
		self.branches = branches 

	def calculate_rates(self, ambientales, loc_irrad, loc_hours_exp, terminal): # this function creates the 
		#random numbers generators to calculate growth rates
		# ambient is a vector with the value of ambient condition at the current time 
		temperature = ambientales[0]
		anomaly = ambientales[1]
		irradiance = loc_irrad # irradiance in kwatts m-2 
		hours_exposition = loc_hours_exp # hours of air exposure at -50 cm

		if loc_hours_exp > 60: # if air exposure exceeds 60 hours the plant wont grow
			def number_of_new_phytomers():
				return 0 

			def length_of_new_phytomers():
				return 0
		else:

			# in python gamma distribution  is parametrized with  a=loc y shape= rate=1/beta
			def number_of_new_phytomers():
				mean_y = 2.154482 # lambda parameter is constant 
				new_phytomers = poisson.rvs(mean_y)
				return new_phytomers

			def length_of_new_phytomers():
				desv = 6.24029427 +(-0.79728514*(anomaly**2))  # the standard  deviation of the distribution
				
				mean_y = (0.47810683 * temperature) + (-0.01180688*(temperature**2)) + (-0.74016238 * (anomaly**2)
					) + (0.15130575 * irradiance) + (-0.13024514 * hours_exposition) #the mean of the distribution

				#using the mean and sd of the distribution calculate the gamma distribution parameters
				length_of_new_phyt = gamma.rvs(a=(mean_y**2)/(desv**2), scale=(1/((mean_y)/(desv**2))), size=1)[0]

				if not terminal: #if the branch is not terminal the internode size is reduced
					length_of_new_phyt = length_of_new_phyt * norm.rvs(loc=0.5642817,scale= 0.08255468)
				return length_of_new_phyt

		return (number_of_new_phytomers, length_of_new_phytomers)
	
	def develop(self, conditions_t, grid_x, grid_y):
		#develop is the method that simulates the growth of a an individual each time step
		#the scheduling  of the develop function is as follows:
		#remove lateral branches with a certain probability
		#add new phytomers to each branch of individual
		#age phytomers and simulate age related mortality
		if len(self.branches)>1:
			lateral_branches = self.branches_here_on(self.branches[0], self.branches[0].phytomers[0], True)#.sort()#lista de los indices de lasr ramas
			lateral_branches = [item for sublist in lateral_branches for item in sublist]
			lateral_branches = sorted(set(lateral_branches),reverse=True) # using set to remove duplicates
			for posible_delete_bran in lateral_branches: #all lateral branches can be separated
				pmuere =  beta.rvs(a=1.672517, b=10.34378,size=1)[0]
				if bernoulli.rvs(pmuere, size=1)[0] == 1: #if there is a succes in separating branch
					index_branch_pointer = self.branches.index(self.branches[posible_delete_bran].origin[0])
					index_phtyo_pointer = self.branches[index_branch_pointer].phytomers.index(self.branches[posible_delete_bran].origin[1])
					self.branches[index_branch_pointer].phytomers[index_phtyo_pointer].delete_branch_out_here()
					self.branch_delete(posible_delete_bran)# separates a branch, the fuctionuses as arguments the pointer 
					#of the branch to remove
		for bran in reversed(self.branches):  # simulating growth
			if bran.active:  # if the meristem of the branch is active
				#calculate local irradiance  and air exposition hours at me meristem grid
				loc_irrad = local_var_from_map(conditions_t[0], grid_x, grid_y,  bran.phytomers[-1])
				loc_hours_exp = local_var_from_map(conditions_t[1], grid_x, grid_y, bran.phytomers[-1])

				(number_of_new_phytomers, length_of_new_phytomers) = self.calculate_rates(conditions_t[2],
				 loc_irrad, loc_hours_exp, bran.terminal())
				new_phytomers= number_of_new_phytomers() #number of new phytomers
				for new in range(new_phytomers): # before adding a phytomer, add lateral branches with a probability
					if bran in self.branches: 
						if bran.terminal() is True:# if the branch is terminal
							pnace =  beta.rvs(a=65.09417, b=127.9993,size=1)[0]#fixed probability
						else: # if the branch is not terminal the probability is reduced
							pnace = beta.rvs(a=65.09417, b=127.9993,size=1)[0] * norm.rvs(loc=0.5642817,scale= 0.08255468)
						if bernoulli.rvs(pnace, size=1)[0] == 1: #if succes in having a branching event
							self.add_branch(length_of_new_phytomers(), self.branches.index(bran)) 	
						#independently to the branching process, add a phytomer
						self.branches[self.branches.index(bran)].add_phyto(length_of_new_phytomers()) 
						self.older_this_branch(1,bran) # ageing of the whole branch

	def refresh(self):  # trough the simulation if a plant turns too small to be alive, or if there are orphaned
	#pointers (pointers that point to instances that are no longer part of an individual) this function deletes
	#these the plants or parts of plants
		for bran in self.branches: #remove pointer pointing at nothing
			for phy in bran.phytomers: 
				if phy.branch_pointer:  # if there is a lateral branch
					if phy.branch_pointer  not in self.branches: #if the pointer points at something not in the indidual
						branch_index = self.branches.index(bran)
						phyt_index = self.branches[branch_index].phytomers.index(phy)
						self.branches[branch_index].phytomers[phyt_index].delete_branch_out_here() # remove 

			if len(bran.phytomers) < 1:
				if (bran.origin[0] in self.branches) and (bran.origin [1] in bran.origin[0].phytomers):
				#if a branch is too small is deleted, and the pointers are also deleted
					index_branch_pointer = self.branches.index(bran.origin[0])
					index_phtyo_pointer = self.branches[index_branch_pointer].phytomers.index(bran.origin[1])
					self.branches[index_branch_pointer].phytomers[index_phtyo_pointer].delete_branch_out_here()
				self.branch_delete(self.branches.index(b))
		

		 #if the individual  has no branches, or if the individual  is formed of a single branch with less phytomers than the minimum
		if (len(self.branches) < 1) or  ((len(self.branches) == 1) and (len(self.branches[0].phytomers) < 5)) :
			if self in meadow:
	
				del meadow[meadow.index(self)]
			if len(meadow) < 1:
				print("The meadow is dead")
	
	def add_branch(self, length, bran):  # bran is the index of the branch instance in the branch list of the individual
		if length > 0 and len(self.branches[bran].phytomers) > 5:  #if the branch is large enough
			last = self.branches[bran].phytomers[-1]  # to get properties of last phytomer of branch
			if self.branches[bran].angle is True:  # the direction of the new branch is determined by the orient variable
				angle = last.orient - 1.169
			elif self.branches[bran].angle is False:
				angle = last.orient + 1.169  
			#create the single new phytomer of the new branch
			new_phyt = phytomer(1, length, last.coord[1], angle)
			self.branches.append(branch([new_phyt],origin=[self.branches[bran],self.branches[bran].phytomers[-1]]))# create
			#the new branch
			self.branches[bran].phytomers[-1].update_branch_out_here(self.branches[-1]) # add the pointer of the new branch to the phytomer of origin
			self.branches[bran].angle_change()  # change the direction of branching 

	def older_this_branch(self, time, bran):  # age all the phytomers in a branch
		for phyt in bran.phytomers:
			phyt.older(time)
		#in a new loop (because eliminating a phytomer would mess the last loop) simulate the probability
		#of dying by age of the first phytomer of a branch (the oldest)
		die_at_age_prob = min(1,max(0,-0.12249192+(0.03649965*bran.phytomers[0].age)))#lineal regression giving the 
		#probability of dying at a certain age, the result may be a number higher than 1 or lower than 0
		#but a probability need to be from 0 to 1, the min and max sets the limits
		if bernoulli.rvs(die_at_age_prob, size=1)[0] == 1:  #if there is a success  at dying
			branch_index = self.branches.index(bran)
			self.phytomer_delete(self.branches[branch_index], self.branches[branch_index].phytomers[0])  

		#the following is a function that returns the branches posterior to a phytomer which is inside a branch
	def branches_here_on(self, bran, phytomer, recursive): 
		#the function returns a list of indexes of branches
		here_on = []  # this is the list that will be filled  
		for phy in bran.phytomers[bran.phytomers.index(phytomer):]:
			if phy.branch_pointer:  # if there is a lateral branch coming out of this phytomer
				branch_index = self.branches.index(phy.branch_pointer)
				here_on.append([branch_index])# add the index 
				if recursive:
					if len(self.branches[branch_index].phytomers) > 0:  # if the branch has at least one phytomer
						here_on = here_on + self.branches_here_on(self.branches[branch_index],
						self.branches[branch_index].phytomers[0],True)

		def filter_empty_list(element):  # function to filter the list from empty elements
			if(element == []):
				return False
			else:
				return True
		here_on = list(filter(filter_empty_list, here_on))  #remove empty lists
		return here_on

	def separate(self, bran, phyt):  # bran and phyt are objects not indexes
		# this function separates the parts of an individual posterior to  breaking point :a phytomer ''phyt'' which 
		# is in a branch ''bran''. the separated individual may become a new individual, and then the separated parts 
		#are deleted from the initial individual
		bran_index = self.branches.index(bran)#the index of the branch of the breaking point
		phyt_index = self.branches[bran_index].phytomers.index(phyt) #the index of the phytomer of the breaking point
		branches_for_new = []# list for the branches that will separate
		pending_pointers = []# list used for updating the pointers after a new individual is created
		list_bye_branch_index = self.branches_here_on(bran,phyt,True)#list of the indexes of the branches posterior to the breaking point
		if len(list_bye_branch_index) > 0: #if there is at least one lateral branch after the breaking point
			#the new individual and the deleted branches of the original will be managed with the index of the branches of the
			#current individual and the index they will have in the new individual
			list_bye_branch_index.insert(0,[self.branches.index(bran)]) # the list will have the index of the terminal branch of
			#the new individual
			list_bye_branch_index = [item for sublist in list_bye_branch_index for item in sublist] #flattens the list of the lateral
			#branches indexes
			list_bye_branch_index = sorted(list(set(list_bye_branch_index)))#set to remove duplicates and then sort
			#if the branch is being removed complete, (the phytomer index is 0) remove the pointer of the phytomer that pointed
			#to this branch 
			if phyt_index == 0:
				self.branches[self.branches.index(bran.origin[0])].phytomers[self.branches[self.branches.index(bran.origin[0])].phytomers.index(bran.origin[1])].delete_branch_out_here()
			for bye_branch in reversed(list_bye_branch_index[1:]): # reversed to not change the following  index numbers when a branch is removed
				#the loop begins at element 1, because element 0 is the terminal branch, which will be processed later
				branches_for_new.insert(0,self.branches[bye_branch]) #add the branches posterior to the breaking point
				pending_pointers.insert(0,
					[list_bye_branch_index.index(self.branches.index(self.branches[bye_branch].origin[0])),
					self.branches[self.branches.index(self.branches[bye_branch].origin[0])].phytomers.index(self.branches[bye_branch].origin[1]),
					list_bye_branch_index.index(bye_branch)])
					#the pointers of the new individual will be updated to match the index in the new individual.
					#the elements of the list pending pointers is a list of lists that contain the index of the branches and pointers
					#that need to be updated: the index of a branch that contains a phytomer with lateral branches, the index of the phytomer
					#the previous branch with a lateral branch, and the index of the branch to which the pointer of the previous phytomer points
				del self.branches[bye_branch] # after extracting the information, delete the branch from the original individual
		#create a new branch that will be the terminal branch of the new individual and add ii to the list 
		branches_for_new.insert(0, branch(bran.phytomers[bran.phytomers.index(phyt):],origin=[]))# create the branch with 
		#all the subsequent phytomers of the breaking  point
		#delete the phytomers or the whole branch from the original individual
		if phyt_index == 0: #if separating the whole branch
			del self.branches[bran_index]
		else: #if only a part of the branch is separating
			del self.branches[bran_index].phytomers[phyt_index:]
		new = zostera(branches_for_new) #create a new individual  
		#the pointers of the new individual need to be updated to point to the branches of the new individual
		#that is, update the points of origin of the branches
		for update in pending_pointers: 
			#each element update has [index of branch of origin, index of phytomer of origin, index of branch coming out of the previous phytomer]
			#the following changes the pointers to the new branches. if the the phytomer of the breakpoint is 0, the indexes of the 
			#phytomer will be the same, else the index of the phytomer will be the original  minus the index of the breakpoint phytomer
			if phyt_index > 0:
				new.branches[update[2]].update_origin([new.branches[update[0]], new.branches[update[0]].phytomers[update[1]-phyt_index]])
				new.branches[update[0]].phytomers[update[1]-phyt_index].update_branch_out_here(new.branches[update[2]])
			else:
				new.branches[update[2]].update_origin([new.branches[update[0]], new.branches[update[0]].phytomers[update[1]]])
				new.branches[update[0]].phytomers[update[1]].update_branch_out_here(new.branches[update[2]])
		self.refresh()
		return new

	def phytomer_delete(self, bran, phy):
		# this function deletes a phytomer of an individual, if that causes to separate two regions, the separate method is called
		#the outcome of this proces depends on the case of the phytomer being deleted:
		bran_index = self.branches.index(bran)
		phyt_index = self.branches[bran_index].phytomers.index(phy) #get pointer
		if not phy.branch_pointer:  # if the pointer is empty 
			if bran.phytomers.index(phy) == 0:
				if bran.terminal():
					# case 1.1: first phytomer in branch (index 0) which is not a branching point and is in the terminal branch
					del self.branches[bran_index].phytomers[phyt_index] #simple delete
				elif bran.terminal() is False:
					# case  1.2: first phytomer in branch (index 0) which is not a branching point and is not in the terminal branch
					if len(bran.phytomers) > 5:#if branch has less 
						new = self.separate(bran, bran.phytomers[bran.phytomers.index(phy) + 1])# separates the rest of the branch
						#and creates a new individual 
						meadow.append(new)#add new individual to meadow list
						try:
							index_branch_pointer = self.branches.index(self.branches[bran_index].origin[0])
							index_phtyo_pointer = self.branches[index_branch_pointer].phytomers.index(self.branches[bran_index].origin[1])
							self.branches[index_branch_pointer].phytomers[index_phtyo_pointer].delete_branch_out_here()
						except ValueError:#the bran fromwhich this phytomer is pointing may have been deleted before, added exception
							pass
					del self.branches[bran_index]
					
			elif bran.phytomers.index(phy) != 0:  # if  phytomer is no the first  in branch
				if phy is bran.phytomers[-1]:
					#case 5, last phytomer in the branch and not a branching point
					del self.branches[bran_index].phytomers[phyt_index]
					self.branches[bran_index].deactivate()# the meristem is gone, so the branch growth is deactivated
				else:
					# caso 2: phytomer is not the first nor last, and is not a branching point
					new = self.separate(bran, bran.phytomers[bran.phytomers.index(phy) + 1])
					meadow.append(new)
					del self.branches[bran_index].phytomers[phyt_index]				
					self.branches[bran_index].deactivate()

		else:  # if phytomer is a branching point
			if bran.phytomers.index(phy) != 0:  # if  phytomer is not the first of the branch
				if phy is bran.phytomers[-1]:
					#case 6, if it the last phytomer in branch and it is a branching point 
					indep_bran_index = self.branches.index(phy.branch_pointer)  
					new = self.separate(self.branches[indep_bran_index], self.branches[indep_bran_index].phytomers[0])
					meadow.append(new)
					del self.branches[bran_index].phytomers[phyt_index]
					self.branches[bran_index].deactivate()

				else:
					#case 3, phytomer that it is neither first nor last in branch and it is a branching point
					indep_bran_index = self.branches.index(phy.branch_pointer)
					new = self.separate(self.branches[indep_bran_index], self.branches[indep_bran_index].phytomers[0])
					meadow.append(new)
					new2 = self.separate(bran, bran.phytomers[bran.phytomers.index(phy) + 1])  
					meadow.append(new2)
					del self.branches[bran_index].phytomers[phyt_index]	
					self.branches[bran_index].deactivate()

			else:# 
				#case 4, phytomer is first in the branch and it is a branching point
				indep_bran_index = self.branches.index(phy.branch_pointer)
				new = self.separate(self.branches[indep_bran_index], self.branches[indep_bran_index].phytomers[0])
				meadow.append(new)
				del self.branches[bran_index].phytomers[phyt_index]				

		self.refresh()


	def branch_delete(self, bran):  # deletes a complete branch, if the branch is large enough it becomes a new individual
		global meadow
		if len(self.branches[bran].phytomers) > 5: #if it has the minimum size for independence
			new = self.separate(self.branches[bran], self.branches[bran].phytomers[0]) #separates and deletes the branch
			# the separate function automatically deletes the parts that are no longer part of the individual
			meadow.append(new)
		else: #if the branch does not have the minumum size it is deleted
			#before it is  deleted, check if any of its phytomers is a branching point, if it is, delete it as well
			if(any([phy.branch_pointer for phy in self.branches[bran].phytomers])):#if it has a phytomer which is a branching point
				to_delete = self.branches_here_on(self.branches[bran],self.branches[bran].phytomers[0],True)
				to_delete.append(bran)
				for bye in (to_delete.sorted(reverse=True)):
					del self.branches[bye]
			else:
				del self.branches[bran]
		self.refresh()

# phytomer is a list of class phytomer instances
class branch:
	def __init__(self, phytomers,origin=[]): 
		self.phytomers = phytomers 
		self.angle = True  # if a new branch will come from the left or right side of the branch
		self.active = True  # indicates if the meristem is active, and thus the branch is growing
		self.origin = origin # list of pointer of the branch from which this branch comes from, and the phytomer
		# from which it developed. if is empty, the branch is terminal
	def __del__(self):
		return 0

	def terminal(self):  #a function to check if the bran is terminal or not
		if self.origin == []:
			is_terminal = True
		else:
			is_terminal = False
		return is_terminal

	def update_origin(self, new_origin=[]): #this function alters the origin variable of a branch,
		#this can be used to indicate that the bran is now terminal (if it has no origin) or when there is 
		#a call to separate and a new individual is created, so the pointer must be updated
		self.origin = new_origin

	def add_phyto(self, length):  # adds a phytomer to a branch of an individual
		#the lenght indicates the internode length of the new phytomer
		if length > 1 and self.active:
			last = self.phytomers[-1]  # to extract information about the last phytomer of the branch
			self.phytomers.append(phytomer(0, length, last.coord[1], last.orient))  # add the new phytomer

	def deactivate(self):
		self.active = False #de activates the growth of a branch, ie, the meristem is no longer alive or active

	def angle_change(self):
		self.angle = not self.angle  # change the branching direction of the next branch that will come out of this branch

	def branch_age(self):
		return len(self.phytomers) #returns  the  number of phytomers in the branch

# phytomer class:
class phytomer:  # the phyomer class represents a phytomer
	def __init__(self, age, length, firstcoord, orient):
		self.age = int(age)  # the age of the phytomer in plastochrons
		self.orient = orient  # the direction in radians toward the phytomer is pointing
		#the following constructor is designed to build phytomers as a succession of partsin a branch, so the coordinates of the 
		#new phytomer are calculated with the coordinates of the node of the previous  phytomer as  the origin of the 
		#internode, as a vertex of a  right triangle, and the internode length as the hypotenuse, so the coordinate of 
		#the node of the newphytomer is calculated with cosine and sine
		self.coord = [firstcoord, [(firstcoord[0] +
				   math.cos(orient) * length), firstcoord[1] +
				   math.sin(orient) * length]] 
		#in the previous line, a  list of two lists of two numbers is created, the firsts are x,y coordinates of the previous node
		#and the second list has the xy coordinates of the new node
		self.branch_pointer = [] # this indicates if the phytomer is the branching point of some branch, in that case, this empty
		#list will have the pointer of that branch, otherwise is empty
		self.length = length

	def older(self, time):
		self.age = self.age + time  # increase the age in plastochrons in some unit of time

	def update_branch_out_here(self, branch_object): #if a branch comes out of this phytomer, the pointer is updated 
		#using this function
		self.branch_pointer = branch_object

	def delete_branch_out_here(self):  # when a branch is lost, the connection  to this phytomer must be eliminated
		#this function elimintes the pointer
		self.branch_pointer = []

#to compare the result of the simulation to the data base od the mathematical ecology laboratory of the CICECSE
#the data of a sample of individual of the meadow are store in a object of this class, this class is used only to arrange the data
#in a certain  order
class ecolmat_data_base():

	def __init__(self):
		self.nodos = pandas.DataFrame(columns=['Date', 'Rhizome', 
	'Internode number', 'Internode length'])

	def add_row(self, meadow, date):
		# when the meadow has more than 20 individuals,  randomly select 20 to make the sample, otherwise, the
		#complete meadow  is sampled 
		#this has two data bases, one with the number of phytomers and internode length of the terminal branch of each individual
		#sampled. And another data base with the number of branches of those same 20 individuals  
		if (len(meadow) > 20):
			sample = np.random.choice(meadow, 20, replace=False)
		else:
			sample = meadow
		
		# the data columns are filled
		tallo_terminal = 0
		for zos in sample:
			tallo_secundario = 0
			for bran in zos.branches:
				if bran.terminal():  # if the branch is terminal
					#terminal
					tallo_terminal = tallo_terminal + 1
					tallo = 'TT%s' % (tallo_terminal)
				else:  # if is lateral branch
					tallo_secundario = tallo_secundario + 1
					tallo = 'TS%s-%s' % (tallo_terminal, tallo_secundario)

			for phyt in zos.branches[0].phytomers:
				#datafrmae contains:
				#'Date','Rhizome', 'Phytomers in terminal branch', 'Internode length'
				num_nodos = len(zos.branches[0].phytomers)
				rizoma = 'R%s' % (tallo_terminal)
				#row has the date of the sample
				row = pandas.Series([date, rizoma, zos.branches[0].phytomers.index(phyt)+1, phyt.length], index=self.nodos.columns)
				self.nodos = self.nodos.append(row, ignore_index=True)
				num_nodos = num_nodos + 1

	def save_data(self, file_name):# this fucntion saves  the data base as a csv file
		place2 = os.getcwd() + os.sep + 'outputs' + os.sep + 'ecol_mat_nod_' + file_name + '.csv'
		place3 = os.getcwd() + os.sep + 'outputs' + os.sep + 'ecol_mat_hoj_' + file_name + '.csv'
		self.nodos.to_csv(place2, encoding='utf-8', index=False)

	#the ambient conditions must be must be in a cvs file with columns for each
	#variable needed and each row is a sequential biweekly set of the observed 
	#or hipotetical state of each variable
	#the file should be ordered as follows
	#      temperat   anomal      irradiance	    hex	0.5 m	 time stamp (posix date)optional
	#time 1    x        y            u                w            fecha
	#time 2    x        y            u                w            fecha
	#time 3    x        y            u                w            fecha
	#time 4    x        y            u                w            fecha
	#  ...
	# to load the  file to make a simulation  
def load_ambient(file_name):  # read the file and organize it as a np array
    route = open(os.getcwd() + os.sep + 'data' + os.sep + file_name, "rt")
    data_csv = csv.reader(route, delimiter=",")
    data = [a for a in data_csv]
    data_np = np.array(data, dtype=float)
    route.close()
    return data_np

def load_meadow(file_name): #load the meadow of the output of a previous simulation named file_name
	place = open(os.getcwd() + os.sep + 'data' + os.sep + file_name, "rb")
	output = pickle.load(place)
	place.close()
	meadow = output
	return meadow


	#fuction to load a world map comprised of and grid of x coordinates, a grid of y coordinates and a 
	#bathymetry of the world. the world components are saved in a pickable file named file_name in the
	#data directory 
def load_map(file_name):
    place = open(os.getcwd() + os.sep + 'data' + os.sep + file_name, "rb")
    output = pickle.load(place)
    place.close()
    (grid_x, grid_y,depth_x_y) = output
    return (grid_x, grid_y, depth_x_y)

    #function to create a map of the variable state at each cell of a world, a map is created for 
	#irradiance and one for hours of exposition to air.
	#the arguments of the function are the world bathymetry, with the depth in mm (or other unit of
	#distance after adjusting the calculation accordingly), a time series of irradiance al the
	#surface, and the hours of exposition at a reference depth of 0.5 m
def variables_map(world, irrad, hours_e):
    #irrad and hours_e must be of the same size
    irradiance = np.empty([len(world), len(world[0]), len(irrad)],
    dtype=float)#create an empty array of same size as world adding a time dimension
    hours_exp = np.empty([len(world), len(world[0]), len(irrad)],
    dtype=float)#again, create another empty array of same size as world adding a time dimension
    alfa = 0.0004  # attenuation coefcient (mm^-1)
    for t in range(len(irrad)):  # make a map for each time , 
        for x in range(len(world)):
            for y in range(len(world[1])):
                #for the irradiance map
                if world[x, y] > 0: #if the is above sea level
                    irradiance[x, y, t] = irrad[t]
                else:#if the cell is is submerged calculate the attenuation with the lambert-beer equation
                    irradiance[x, y, t] = (irrad[t] * math.exp(abs(world[x, y]/100) * -alfa)) 
                #for the hours of exposition to air in each cell
                if world[x, y] < -566:  #if the sea level is bellow a limit (this one is for the punta banda estuary)
                #the cell is always submerged
                    hours_exp[x, y, t] = 0
                elif world[x, y] > 2110:  #if the cell is higher than this level  it is always exposed to air
                    hours_exp[x, y, t] = 360
                #if the depth is between the previous limits make a lineal interpolation using the two limits and the
                # hours exp at the 0.5m reference  value
                else:
                    hex_for_interp=[0, hours_e[t], 360]  #known  hours of exposition at the following depths
                    prof_for_interp=[-566, 500, 2110] #the depths
                    interpolater = interp1d(prof_for_interp, hex_for_interp, kind='linear') #use linear  interpolation
                    #to 
                    hours_exp[x, y, t] = interpolater(world[x, y])#use the lineal interpolator at the cell 
    var_maps = [irradiance, hours_exp]#put together the generated maps
    return var_maps

def create_zosteras_from_csv(file_name): # this function creates a meadow from the information of a csv file.
	#in the csv file there is information to create a number induvials formed by only one terminal branch.
	#the csv has information in columns  which contain in order: first coordinate of oldest pytomer of a branch,  y coordinate,
	#orientation of the branch in radians, the internode length  of each internode in the branch
	place = os.getcwd() + os.sep + 'data' + os.sep + file_name
	internode_lengths = (pandas.read_csv(place, header=None).values)
	created_zosteras =[]
	# in a loop create individuals
	for individual in range(len(internode_lengths)): 
		fila = internode_lengths[individual] #extract the data
		lista = [phytomer(0, fila[3], [fila[0], fila[1]], fila[2])]  # create the first phytomer
		rama = branch(lista)#create the branch with  one phytomer
		created_zosteras.append(zostera([rama]))  # create the branch and make it terminal
		for internod in range(4,len(fila)):  #for each internode lenght column
			if not np.isnan(fila[internod]):
				for internode_ages in created_zosteras[individual].branches[0].phytomers:
					internode_ages.older(1)	# age the internodes  one plastochron, as they should be ordered from older to youngest
				created_zosteras[individual].branches[0].add_phyto(fila[internod]) # add phytomer
	print("Meadow created with %d individual from the file %s" %( len(internode_lengths), file_name))
	return created_zosteras

def local_var_from_map(var_map, grid_x, grid_y, phytomer):
	# function to calculate the irradiance oh hours of exposition at a certain cell of the world  grid
	loc_var = griddata((grid_y.ravel(), grid_x.ravel()), (var_map[:, :]).ravel(), (phytomer.coord[0][0], phytomer.coord[0][1]),
	method='nearest') 
	return loc_var

# the next part is to plot the meadow in a time t
points = []  # list to save the points to plot (the phytomers)
ages = []  # the list of each phytomer age
lengths = [] #the lenght of the internodes
zos_num = []# the number of individual in the meadow
#date_str = []#the date of t as a character
demogra = [[], [], [], []]  # demographic data to be plotted

meadow = []  # the list of the zostera class instances, represents  the individuals alive

def data_saver(meadow, t, date_str):
	#this fuctions is to fill the previous lists, each time  a new set of measures will be added to the list
	points.append([])
	ages.append([])
	zos_num.append([])
	lengths.append([])#add the space to be filled with the observations at t

	demogra[0].append(len(meadow))  # population size
	demogra[2].append(1)  # empty space, could be used for density of individuals by area
	demogra[3].append(date_str)#the date of t as a character
	#iniciar contador de cuantos tallos hay
	brach_count = 0
	for zos in meadow:
		brach_count = brach_count + (len(zos.branches))
		for bran in zos.branches:
			for phyt in bran.phytomers:
				points[t].append(phyt.coord)
				ages[t].append(phyt.age)
				zos_num[t].append(meadow.index(zos))
				lengths[t].append(phyt.length)
	demogra[1].append(brach_count)
	

def simulation(meadow, ambientales, var_maps, grid_x, grid_y, save, *args):
	#the loop of a simulation, represents the flow of time biweekly.
	#each element of the ambient time series the individuals develop accordingly to the
	#conditions of those two weeks.
	#This functions needs the time series ambient, initial meadow, the map of the variables
	#irradiance and hours of exposition at each cell of the world at each time, the grid of the 
	#save is a boolean specifying if the result of the simulation is to be saved as a file in the mathematical ecology lab format
	#and optionally an extra argument, a character for the name of the output data.

	if save is True: #initialize de datafile to be filled
		data_ecolmat = ecolmat_data_base()

	for t in tqdm(range(len(ambientales))):#tqdm is to show a progress bar of the loop
		conditions_t = [var_maps[0][:,:,t], var_maps[1][:,:,t],ambientales[t]] #joins in a list the map of irradiance an hours of exp
		#at time t, and the rest of the ambient information at time t

		for zos in meadow:
			zos.develop(conditions_t, grid_x, grid_y)  # each individual develops according to the ambient conditions
		#sample and record the state of the sampled individuals at the end of the time, this represents the end of the 
		#sampling interval
		try:
			data_saver(meadow, t, datetime.fromtimestamp(ambientales[t][4]).strftime("%d/%m/%Y"))#  save data for plots
			if save is True:  # to save data for the csv data file 
				data_ecolmat.add_row(meadow, datetime.fromtimestamp(ambientales[t][4]).strftime("%d/%m/%Y"))
		except IndexError:
			data_saver(meadow, t, str(t))
			if save is True: 
				data_ecolmat.add_row(meadow, str(t))

	if save is True: #at the end of the simulation the csv file is saved (if requested)
		#and the outputs are also saved using the pickle module
		if not (os.path.exists(os.getcwd() + os.sep + 'outputs')):
			os.mkdir('outputs')
		file_name = "output.dat"
		if args:
			file_name = args[0]
		place1 = open(os.getcwd() + os.sep + 'outputs' + os.sep + file_name + '.dat',
			 "wb")
		pickle.dump((meadow), place1)  
		place1.close()
		#saves the csv  file
		data_ecolmat.save_data(file_name)

		print(("Outputs saved as %s" % file_name))
	return (points, ages, lengths, demogra, zos_num)

def main(inputs, make_plot, seed = 26): #the main function
	rn.seed(seed)  # set the random seed to use in the random module
	np.random.seed(seed)  #and in the numpy module
	#the inputs list has strings of the name of the files to use in the simulation in order: meadow, ambient, world.
	#this files must be in a 'data' directory in the directory of this file. the initial meadow can be
	#the meadow obatined from a previous simulation, in that case the first argument must be the name
	#of sich file, to create a meadow using a csv file set te argument to be a string that does not match
	#any file in the inptus directory.
	#for the demonstration simulation use the files:
	#for the ambient condtitions: 'ambient_2000.csv'
	#for the csv file of the initial branches: 'initial_branches_2000.csv'
	#for the world: 'cannal_200m_broad_4m_prof.data', this file can be created in 'sample_world.py' 
	#inputs[1:] = ["sample_founding_rhizomes_2000.csv", "sample_environment2000.csv", "cannal_200m_broad_4m_prof.dat"]

	print ("Preparing simulation")
	#inputs are Initial, Environment, World
	#STEP  1: create a meadow of load one
	#load a meadow if there is a pickable file in in the data directory
	if inputs[0].endswith(".csv"):
		for m in create_zosteras_from_csv(inputs[0]):  # set the name of  the csv file here
			meadow.append(m)
		print("Meadow succesfully created")
	elif os.path.isfile(os.getcwd() +  os.sep + 'data' + os.sep + inputs[0]): #if there is a file that matches the meadow file name
		for m in load_meadow(inputs[0]):
			meadow.append(m)
		print("Meadow loaded from file: " + inputs[0])
	else:
		raise NameError("Invalid Initial file")
		
	#STEP 2 import ambiental condition and world
	if os.path.isfile(os.getcwd() +  os.sep + 'data' + os.sep + inputs[1]):
		ambiental = load_ambient(inputs[1])
		print("Ambient conditions time series loaded from file: " + inputs[1])
	else:
		raise NameError("Invalid Environment file")

	#load world map
	if os.path.isfile(os.getcwd() +  os.sep + 'data' + os.sep + inputs[2]):
		(grid_x, grid_y, depth_x_y) = load_map(inputs[2])
		print("World loaded from file: " + inputs[2])
	else:
		raise NameError("Invalid World file")

	#grid_x and grid_y is the mesh grid of the coordinates of the cells of the world, dept_x_y has the depth of those
	#same cells
	#create map of variables at each cell at each time
	var_maps = variables_map(depth_x_y, ambiental[:, 2],ambiental[:, 3])
	# var_maps has two maps, the first of irradiance and the second of hours of exposition
	print("Ambient and world ready")
	
	#STEP 3. make a simulation iterating over the length of the ambient time series
	start_time = time.time()#to print the duration of the simulation 

	print(("Simulating..."))

	output_name = '_'.join([file.split('.')[0] for file in inputs]) + "_" + str(seed)# the name of the output is the combination of the input files

	output = simulation(meadow, ambiental, var_maps, grid_x, grid_y, True, output_name)

	print("Simulaction finished. Time elapsed : %s seconds" % round((time.time() - start_time), 1))
	#Step 4 optionally make an animation of the simulation and or save  it

	plot_meadow(output, var_maps[0], [grid_x, grid_y, depth_x_y], save=make_plot, show=False)# the plotting function
	#the arguments of the plot meadow function are the output of a simulation, a map of the variable to plot over
	#, or the depth, and as a list the coordinates of y grid, coordinates of x grid and the depth, save and show are
	#booleans that indicate if an animation should be made and showed, and if it should be saved as mp4

	print("The final population size is %s individuals" %len(meadow))

	print("~(^w^)~")#succes

if __name__ == "__main__":
	make_plot_arg = sys.argv[4].lower() == 'true'
	if len(sys.argv) == 6:#if the four arugments (plus file name) are specified
		main(inputs = sys.argv[1:4], make_plot = make_plot_arg, seed = int(sys.argv[5]))
	else:
		main(inputs = sys.argv[1:4], make_plot = make_plot_arg)#else use the default seed
