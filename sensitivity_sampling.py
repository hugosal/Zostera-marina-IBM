#o-*- coding: utf-8 -*-
#Hugo Salinas 2020
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
from tqdm import tqdm


#the following simlation is an adaptation of the one in zostera_model.py to be used in the sensitibvity analysis
#this function will be executed once for each of the k rows of the sample spacematrix.
#fila_de_parametros is a set of paramter values, see bellow.
def simulacion_con_fila(fila_de_parametros):
	meadow=[]

	for m in create_zosteras_from_csv('founding_rhizomes_2000.csv'):  # cambiar nombre archivo aqui
		meadow.append(m)

	ambiental = load_ambient("environment_2000.csv")

	(grid_x, grid_y, depth_x_y) = load_map("cannal_200m_broad_4m_prof.dat")

	var_maps = variables_map(depth_x_y, ambiental[:, 2],ambiental[:, 3])

	output = simulation(meadow, ambiental[0:9], var_maps, grid_x, grid_y, fila_de_parametros)
	
	return output

class zostera:
	def __init__(self, branches):  # branches es un lista de ramas
		self.branches = branches  # lista objetos rama

	def calculate_rates(self, ambientales, loc_irrad, loc_hours_exp, terminal, fila_de_parametros):    # para calcular las tasas de crecimiento
		# ambiente es un vector con la temp, anomalia, etc. 
		(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10,
		theta11, theta12, theta13, theta14, theta15, theta16) = fila_de_parametros
		temperature = ambientales[0]
		anomaly = ambientales[1]
		irradiance = loc_irrad
		hours_exposition = loc_hours_exp
		#print("global ")
		#print(ambientales[3])
		#print("local")
		#print(hours_exposition)


		if loc_hours_exp > 60: # el dato mas grande de horas es de 60, y eso significa que no se vale extrapolar a mas
			#alla de ese valor, ademas, si estan tantotiempo expusto no deberian recer, asi quepongo estoasi
			print("estoy en tierra")
			def number_of_new_phytomers():
				return 0 

			def length_of_new_phytomers():
				return 0
		else:
			#los parametros anteriores se calcularon bayesianamente en  regresion_crec_long.r, para la simulacion 
			#se va a generar un numero aleatorio de una distribucion cuya forma esta dad por parametros que siguen
			#la distribucion gamma aquise parametriza con a=loc y shape= rate=1/beta
			def number_of_new_phytomers():
				#generador poisson poisson.rvs(1, size=10)
				mean_y = theta1 # lambda constante 
				new_phytomers = poisson.rvs(mean_y)
				return new_phytomers

			def length_of_new_phytomers():
				#x,x2,y,u,w,uw
				#la desviacion de la distribucion varia con la anomalia cuadrada:
				desv = theta2 +(theta3*(anomaly**2))
				#la media de longitud es:
				mean_y = (theta4 * temperature) + (theta5*(temperature**2)) + (theta6 * (anomaly**2)
					) + (theta7 * irradiance) + (theta8 * hours_exposition)
				#length_of_new_phyt =  max(0, norm.rvs(loc=mean_y,scale=desv, size=1)[0])
				length_of_new_phyt = gamma.rvs(a=(mean_y**2)/(desv**2), scale=(1/((mean_y)/(desv**2))), size=1)[0]
				if not terminal:
					#el siguiente es para reducir el tamanio de los internodos de las ramas laterales, el valor
					#se calculo con una muestra pequenia, los numeros viene de longitud_rizoma_vs__lateral.R
					length_of_new_phyt = length_of_new_phyt * norm.rvs(loc=theta9,scale= theta10)
				#print(length_of_new_phyt)
				return length_of_new_phyt

		return (number_of_new_phytomers, length_of_new_phytomers)
	
	def develop(self, conditions_t, grid_x, grid_y, fila_de_parametros):
		(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10,
		theta11, theta12, theta13, theta14, theta15, theta16) = fila_de_parametros
		#el nuevo schedulling de develop es el siguiente:
		#dentro de un ciclo for para una rama 1. quitar ramas hereon  recursiva con prob fija que 
		#viene de regresion_crec_rama_bayes_ver_2.R, despues agregar con un for los fitomeros que prediga la funcion
		#number_of_new_phytomer, pero cada vez que se agregue uno, antes de que salga puedesalir una ramacon prob fija,
		#ya que se agrega un fitomero auna rama se envejece la rama y con prob en funcion de edad se muere ultimo fitomero
		#de la rama
		if len(self.branches)>1:
			lateral_branches = self.branches_here_on(self.branches[0], self.branches[0].phytomers[0], True)#.sort()#lista de los indices de lasr ramas
			lateral_branches = [item for sublist in lateral_branches for item in sublist]
			lateral_branches = sorted(set(lateral_branches),reverse=True) #set paraque no se repitan
			for posible_delete_bran in lateral_branches: #las ramas pueden morir
				#pmuere = min(max(0,norm.rvs(loc=0.2460341, scale=0.09731517,size=1)),1)#esta es laultima.. peronoda buen resultado
				#la distribucion beta se puede parametrizar conmedia y desv de la distribucion como
				#alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
  				#beta = alpha * (1 / mu - 1)
				#pmuere = 0.2040021#0.4308347
				pmuere =  beta.rvs(a=theta11, b=theta12,size=1)[0]
				if bernoulli.rvs(pmuere, size=1)[0] == 1: #esta prob es fija, si la rama debe morir la borro
					index_branch_pointer = self.branches.index(self.branches[posible_delete_bran].origin[0])
					index_phtyo_pointer = self.branches[index_branch_pointer].phytomers.index(self.branches[posible_delete_bran].origin[1])
					self.branches[index_branch_pointer].phytomers[index_phtyo_pointer].delete_branch_out_here()
					self.branch_delete(posible_delete_bran)# las separa, branch_deleteo toma como argumento indice de rama
					# se elimina el pointer de  la rama, esdecir, se quitar el pointer del fitomero de donde salia la rama
				

		for bran in reversed(self.branches):  # ahora para crecer y ramificar
			if bran.active:  # si el meristemo de la rama esta activo:si esta creciendo 
				#calcular irradiancia y hex local en el punto del fitomero que crece 
				loc_irrad = local_var_from_map(conditions_t[0], grid_x, grid_y,  bran.phytomers[-1])
				loc_hours_exp = local_var_from_map(conditions_t[1], grid_x, grid_y, bran.phytomers[-1])

				(number_of_new_phytomers, length_of_new_phytomers) = self.calculate_rates(conditions_t[2],
				 loc_irrad, loc_hours_exp, bran.terminal(), fila_de_parametros)
				new_phytomers= number_of_new_phytomers() # cuantos van a salir
				for new in range(new_phytomers): #antes de agregar los fitomeros que saldrian, puede salir una rama de este 
					if bran in self.branches:  # puede que la rama ya no este por algun otro metodo que la borra
						if bran.terminal() is True:#si es terminal
							pnace =  beta.rvs(a=theta13, b=theta14,size=1)[0]#0.3371122
							#pnace = norm.rvs(loc=0.3357790, scale=0.03192818, size=1)#esta es laultima,perono da buen resultado
						else: 
							pnace = beta.rvs(a=theta13, b=theta14,size=1)[0] * norm.rvs(loc=theta9,scale= theta10)
						if bernoulli.rvs(pnace, size=1)[0] == 1: # la rama se ramifica con probabilidad fija
							self.add_branch(length_of_new_phytomers(), self.branches.index(bran)) #sale rama	
						#independiente de si se ramifica o no se agregafitomero 
						self.branches[self.branches.index(bran)].add_phyto(length_of_new_phytomers()) #se agregafitomero arama central
						self.older_this_branch(1,bran, fila_de_parametros) #se envejece ramainicial

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

	def older_this_branch(self, time, bran, fila_de_parametros):  # age all the phytomers in a branch
		(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10,
		theta11, theta12, theta13, theta14, theta15, theta16) = fila_de_parametros
		for phyt in bran.phytomers:
			phyt.older(time)
		#in a new loop (because eliminating a phytomer would mess the last loop) simulate the probability
		#of dying by age of the first phytomer of a branch (the oldest)
		die_at_age_prob = min(1,max(0,theta15+(theta16*bran.phytomers[0].age)))
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

def load_ambient(file_name):  # read the file and organize it as a np array
    route = open(os.getcwd() + os.sep + 'data' + os.sep + file_name, "rt")
    data_csv = csv.reader(route, delimiter=",")
    data = [a for a in data_csv]
    data_np = np.array(data, dtype=float)
    route.close()
    return data_np


def load_meadow(file_name): #carga las zosteras de un output de una simulacion anterior
	place = open(os.getcwd() + os.sep + 'data' + os.sep + file_name, "rb")
	output = pickle.load(place)
	place.close()
	meadow = output
	return meadow

def create_zosteras_from_csv(file_name): # crea zosteras con internododos detamanio x espcificados en el archivo
	#el csv debe tener los datos de un rizoma por fila, las columnas son las primeras 2 las coordenadas, la 
	#tercera la orientacion, y de la 4al final las longitudes de los internodos
	#cargarel archivo
	place = os.getcwd() + os.sep + 'data' + os.sep + file_name
	internode_lengths = (pandas.read_csv(place, header=None).values)
	created_zosteras =[]
	#phytomer con  rhizome, branch, prev_fit, age, length, firstcoord, orient):
	for individual in range(len(internode_lengths)):  #crear una zostera por cada fila
		fila = internode_lengths[individual] #saco los datos de la long de cada internodo
		lista = [phytomer(0, fila[3], [fila[0], fila[1]], fila[2])]  # creo el primer fitomero
		rama = branch(lista)
		created_zosteras.append(zostera([rama]))  #creo  la rama, la hago terminal , la hago zostera y la 
		#meto a lista en una sola linea, ahora a esa zostera le voy a a gregar los fitomeros que faltan
		for internod in range(4,len(fila)):  #por  cada fitomero que falta agregar
			if not np.isnan(fila[internod]):
				for internode_ages in created_zosteras[individual].branches[0].phytomers:
					internode_ages.older(1)		#quiero que tengan la edad adecuada, osea en escalera, el mas viejo al final
				created_zosteras[individual].branches[0].add_phyto(fila[internod]) #  a la zostera en cuestion le agrego el fitomero
	#print("Pradera creada con %d individuos del archivo %s" %( len(internode_lengths), file_name))
	return created_zosteras

def local_var_from_map(var_map, grid_x, grid_y, phytomer):
	# funcion para calcular la irradiancia o hex en un punto especifico, dado que
	#el fitomero esta en un cierto punto de la matriz
	loc_var = griddata((grid_y.ravel(), grid_x.ravel()), (var_map[:, :]).ravel(), (phytomer.coord[0][0], phytomer.coord[0][1]),
	method='nearest') #si, le estoy dando primero la grilla de y y luego la dey,porque asi esta funcionando... no se por
	return loc_var

def data_saver(meadow ):
	#crear los espacio para llenarse con los datos de cada tiempo, meadow es meddow, t es int con el indice de
	#las coniciones ambientales del momento, date_str es el nombre de la fecha en fomrato "d/m/Y"
	#hay dos tipos de variables a guardar: las que se guardan una vez por tiempo (numero
	#de zosteras, temp, etc), y las que se guardan por cada fitomero (long, edad, etc: las 
	#que se guardan por fitomero:
	lengh_list = []
	for zos in meadow:
		this_ind=[]
		for phyt in zos.branches[0].phytomers:
			this_ind.append(phyt.length)
		lengh_list.append(sum(this_ind))
	return sum(lengh_list)/len(lengh_list)
	
def simulation(meadow, ambientales, var_maps, grid_x, grid_y, fila_de_parametros):
	#simlacion es el ciclo principal que itera sobre el tiempo. cada tiempo t de la 
	#serie de tiempo de concidionesambientales se desarrollan las zosteras de acuerdo a ese
	#tiempo.
	# simulation necesia una pradera, una serie de tiempo de ambiente, el mundo
	#la malla en x y y, la indicacion se si  guardar o no, y opcinal el nombre del
	#archivo si se va a guardar (args)
	#ambientales vienen en dimensiones (t, 5), en orden esas 5 variables son: temp, anom, irrad, hex y timestamp
	#var_maps  es una lista con dos np array de (x,y,t) donde x y y son las dimensiones espaciales y t la temporal, el primer
	#elemento de la lista varmaps es el array con los datos de irrad en (x,y) y la seghunda las hexen (x,y) 
	

	for t in range(len(ambientales)):#tqdm(range(len(ambientales))): 
		#plot_zosters(meadow) #grafica de como empezo
		#print(len(meadow)) #activar alguna de estas quita el chiste a tqdm
		
		conditions_t = [var_maps[0][:,:,t], var_maps[1][:,:,t],ambientales[t]] #esta lista tiene los mapas de ambientales 
		#al tiempo t, y las condiciones al tiempo t, se usa la misma para todas las zosteras:todas estan expuestas a lo mismo
		
		for zos in meadow:
			zos.develop(conditions_t, grid_x, grid_y, fila_de_parametros)  # el develop necesita la lista de condiciones altiempo y las mallas
			#en x y y parapoder interpolar.

		#se guarda la pradera como estaba la pradera al terminar el desarrollo, osea, como si fuera
		#al momento de la recaptura de los rizomas
		mean_len = data_saver(meadow)  # de ley para hacer grafica
	return mean_len

#the sample space of the parameters was obtained with latin hypercube method.
#This space corresponds to entries ina matrix of n rows and k columns. Each n row is a set of paramter value,
#each of the k columns are the values of a parameter
#the range of the parameter values to be tested is plus/minus 30 percent of the actual parameter value

rn.seed(26)
np.random.seed(26)

#to fix the size of the matrix with n rows and k columns
num_param = 16 #there are 16 parameters to be studied in the model
extra = 84 #Marino2009 suggest a minimum of n to be num_param+1, but is it better to add more
n_row =  num_param + extra # the final number of n rows

#to get the value of the entries of the matrix.
rangos = np.zeros((num_param, 4), dtype=float ) #aqui va a estar en columnas el maximo de cada 
#rango, y el valor original, el rango es el valor original mas menos 0.3 de ese valor
rangos[:,1] = [2.15, 6.24, -0.79, 0.4781, -0.011, -0.740, 0.1513, -0.1302, 0.5642817,
	 0.08255468,1.672517,10.34378, 65.09417,127.9993, -0.12249192, 0.03649965]#set of the actual parameter values
rangos[:,0]=rangos[:,1]-(rangos[:,1]*0.3) # inf limit
rangos[:,2]=rangos[:,1]+(rangos[:,1]*0.3) # sup limit
rangos[:,3]=(rangos[:,2]-rangos[:,0])/n_row # jump size

parametros = np.zeros((n_row, num_param), dtype=float ) #inititalize the matrix

for f in range(n_row): 
	for p in range(num_param):
		parametros[f, p] = rn.uniform((rangos[p,0]+(f*rangos[p,3])), (rangos[p,0]+((f+1)*rangos[p,3]))) # set a value in the range
for p in range(num_param):#shuffle the columns
	np.random.shuffle(parametros[:,p])
np.savetxt('latin_hipercube.csv', parametros, delimiter=',') #save file, this file will be used in the computation of the correlation

#simulation using each of rows of parameter valueas input
resultado = np.zeros((num_param+extra,1))
for fila in tqdm(range(len(parametros[:,1]))):
	rn.seed(12)# set the same seed in all simulations
	np.random.seed(12) 
	lengths = []
	meadow = [] 
	resultado[fila] = simulacion_con_fila(parametros[fila,:])

np.savetxt('outputs_for_latin.csv', resultado, delimiter=',')#save the  output meassures in a file