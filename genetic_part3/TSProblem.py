# -*- coding: utf-8 -*-
"""
Created on Feb 18 2022
Last edit on April 02 2022

@author: tdrumond & agademer & ndenier & hgu

Template file for your Exercise 3 submission 
(GA solving TSP example)
"""
from GASolver_module import GAProblem
import cities_module as cm
import random

class TSProblem(GAProblem):
    """Implementation of GAProblem for the traveling salesperson problem"""
    def __init__(self, possible_genes, threshold_fitness):
        super().__init__(possible_genes, threshold_fitness)
        self.nbCities = cm.countCities("cities.txt")

    def generate_chromosome(self):
        chromosome = self._possible_genes
        random.shuffle(chromosome)
        
        return chromosome

    def fitness_function(self, chromosome):
        return -cm.roadLength(cities, chromosome) 
        
    def crossing_function(self, population, selection_rate): # return new chromosome
        # Select parent 1 randomly
        parent1 = population[0] # Default value in case no parent is selected
        for i in population:
            number = random.random()
            if (number < selection_rate):
                parent1 = i.chromosome
                parent1_index = population.index(i)
                break
        # Select parent 2 randomly andd different from parent 1
        parent2 = population[1] # Default value in case no parent is selected
        for i in population:
            number = random.random()
            if (number < selection_rate and population.index(i) != parent1_index):
                parent2 = i.chromosome
                break
     
        # Create new individual from 2 parents
        cut = round(self.nbCities/2)
        new_chromosome = parent1[:cut]
        for c in parent2[cut:]:
            if c not in new_chromosome:
                new_chromosome.append(c)
        for c in self._possible_genes:
            if c not in new_chromosome:
                new_chromosome.append(c)

        return new_chromosome

    def mutation_function(self, population, mutation_rate):
        # Mutation for _mutation_rate of the population
        for i in population:
            number = random.random()
            if (number < mutation_rate):
                # A random value is assigned to a random gene of the chromosome
                city1 = random.randint(0,self.nbCities-1)
                city2 = random.randint(0,self.nbCities-1)
                while (city1 == city2):
                    city2 = random.randint(0,self.nbCities-1)
                # swap two cities
                i.chromosome[city1], i.chromosome[city2] = i.chromosome[city2], i.chromosome[city1]
                i.fitness = -cm.roadLength(cities, i.chromosome) # The fitness is recomputed for the mutated chromosome

        return 0


if __name__ == '__main__':

    from GASolver_module import GASolver

    cities = cm.loadCities("cities.txt")
    problem = TSProblem(possible_genes=cm.defaultRoad(cities), threshold_fitness=-350)
    solver = GASolver(problem)
    solver.resetPopulation()
    solver.evolveUntil()
    best = solver.getBestIndividual()
    
    print(best.chromosome, best.fitness)
    #cm.drawCities(cities, solver.getBestIndiv().chromosome)
