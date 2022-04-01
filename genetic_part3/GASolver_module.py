# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 2022

@author: tdrumond & agademer & ndenier & hgu

Template file for your Exercise 3 submission 
(generic genetic algorithm module)
"""
import random


class Individual:
    """Represents an Individual for a genetic algorithm"""

    def __init__(self, chromosome: list, fitness: float):
        """Initializes an Individual for a genetic algorithm 

        Args:
            chromosome (list[]): a list representing the individual's chromosome
            fitness (float): the individual's fitness (the higher, the better the fitness)
        """
        self.chromosome = chromosome
        self.fitness = fitness

    def __lt__(self, other):
        """Implementation of the less_than comparator operator"""
        return self.fitness < other.fitness

    def __repr__(self):
        """Representation of the object for print calls"""
        return f'Indiv({self.fitness:.1f},{self.chromosome})'


class GAProblem:
    """Defines a Genetic algorithm problem to be solved by GASolver"""
    def __init__(self, possible_genes: list, threshold_fitness: int):
        """
        """
        self._possible_genes = possible_genes
        self._threshold_fitness = threshold_fitness

        def generate_chromosome(self): # return chromosome
            return 0

        def fitness_function(self, chromosome):
            return 0
        
        def crossing_function(self): # return new chromosome
            return 0

        def mutation_function(self):
            return 0


class GASolver:
    def __init__(self, problem: GAProblem, selection_rate=0.5, mutation_rate=0.1):
        """Initializes an instance of a GASolver for a given GAProblem

        Args:
            problem (GAProblem): GAProblem to be solved by this GASolver
            selection_rate (float, optional): Selection rate between 0 and 1.0. Defaults to 0.5.
            mutation_rate (float, optional): mutation_rate between 0 and 1.0. Defaults to 0.1.
        """
        #self._possible_genes = problem._possible_genes
        self._problem = problem
        self._selection_rate = selection_rate
        self._mutation_rate = mutation_rate
        self._population = []
        self.generation = 0

    def resetPopulation(self, pop_size=50):
        """ Initialize the population with pop_size random Individuals 

        Args:
            pop_size (int, optional): size of the population (positive, minimum related to the selection_rate). Default to 50.
        """        
        self._population = [] # Empty the population

        # Generate pop_size random chromosomes
        for i in range(pop_size):
            chromosome = self._problem.generate_chromosome()
            #chromosome = self._problem._possible_genes
            #random.shuffle(chromosome)
            
            fitness= self._problem.fitness_function(chromosome) # Compute the fitness of the chromosome

            new_individual = Individual(chromosome=chromosome, fitness=fitness) # Create new individual
            self._population.append(new_individual) # Add it to the population


    def evolveForOneGeneration(self):
        """ Apply the process for one generation : 
            -	Sort the population (Descending order)
            -	Remove x% of population (less adapted)
            -   Recreate the same quantity by crossing the surviving ones 
            -	For each new Individual, mutate with probability mutation_rate 
                i.e., mutate it if a random value is below mutation_rate"""
        
        # Sort the population (descending order)
        self._population.sort(reverse=True)

        # Remove _selection_rate of the population
        number_individuals = len(self._population)
        number_kept = round(self._selection_rate*number_individuals)
        self._population = self._population[:number_kept]

        #Crossing (select 2 parents then create a new individual)
        for n in range(number_individuals-number_kept):
            new_chromosome=self._problem.crossing_function(population=self._population, selection_rate=self._selection_rate) 
            new_individual = Individual(chromosome=new_chromosome, fitness=self._problem.fitness_function(new_chromosome))
            self._population.append(new_individual)

        # Mutation
        self._problem.mutation_function(self._population, self._mutation_rate)

        # Increment the generation count
        self.generation+=1 


    def showGenerationSummary(self):    
        """ Print some debug information on the current state of the population """
        print('Generation: ',self.generation)
        print('chromosome\tfitness')
        for i in self._population:
            print(i.chromosome,'\t',i.fitness)

    def getBestIndividual(self):
        """ Return the best Individual of the population """
        self._population.sort(reverse=True)
        best=self._population[0]
        return best

    def evolveUntil(self, max_nb_of_generations=500):
        """ Launch the evolveForOneGeneration function until one of the two condition is achieved : 
            - Max nb of generation is achieved
            - The fitness of the best Individual is greater than or equal to
              threshold_fitness
        """
        fitness_max = self.getBestIndividual().fitness
        while (self.generation < max_nb_of_generations and fitness_max < self._problem._threshold_fitness):
            fitness_max = self.getBestIndividual().fitness
            self.evolveForOneGeneration()