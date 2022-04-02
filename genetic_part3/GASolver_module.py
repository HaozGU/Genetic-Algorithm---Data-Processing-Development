# -*- coding: utf-8 -*-
"""
Created on Feb 18 2022
Last edit on April 02 2022

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
        """Initializes a problem

        Args: 
            possible_genes (list): a list of all possible genes
            treshold_fitness (int): a satisfying value of fitness to stop the algorithm
        """
        self._possible_genes = possible_genes
        self._threshold_fitness = threshold_fitness

        # functions different for each problem:
        def generate_chromosome(self):
            """ Generate a new chromosome
                Return value: 
                    chromosome (list): a random chromosome
            """
            return 0

        def fitness_function(self, chromosome):
            """ Calculate the fitness value for a given chromosome
                Args:
                    chromosome (list)
                Return value:
                    fitness (int): the fitness value of the chromosome
            """
            return 0
        
        def crossing_function(self):
            """ Select two parents and then generate a new chromosome by merging the parent's chromosomes
                Return value: 
                    new_chromosome (list): a new chromosome generated from the chromosomes of two different individuals (parents)
            """
            return 0

        def mutation_function(self):
            """ Mutate randomly some chromosomes for some randomly chosen individuals of a population
            """
            return 0


class GASolver:
    def __init__(self, problem: GAProblem, selection_rate=0.5, mutation_rate=0.2):
        """Initializes an instance of a GASolver for a given GAProblem

        Args:
            problem (GAProblem): GAProblem to be solved by this GASolver
            selection_rate (float, optional): Selection rate between 0 and 1.0. Defaults to 0.5.
            mutation_rate (float, optional): mutation_rate between 0 and 1.0. Defaults to 0.2.
        """
        self._problem = problem
        self._selection_rate = selection_rate
        self._mutation_rate = mutation_rate
        self._population = [] # list of all the individuals
        self.generation = 0 # to count the generations

    def resetPopulation(self, pop_size=50):
        """ Initialize the population with pop_size random Individuals 

        Args:
            pop_size (int, optional): size of the population (positive, minimum related to the selection_rate). Default to 50.
        """        
        self._population = [] # Empty the population

        # Generate pop_size random chromosomes:
        for i in range(pop_size):
            # Generate a random chromosome
            chromosome = self._problem.generate_chromosome()
            # Compute the fitness value of the chromosome
            fitness= self._problem.fitness_function(chromosome) 
            # Create a new individual with the chromosome and the associated fitness value
            new_individual = Individual(chromosome=chromosome, fitness=fitness) 
            # Add it to the population
            self._population.append(new_individual) 


    def evolveForOneGeneration(self):
        """ Apply the process for one generation : 
            -	Sort the population (Descending order)
            -	Remove x% of population (less adapted)
            -   Recreate the same quantity by crossing the surviving ones 
            -	For each new Individual, mutate with probability mutation_rate 
        """
        
        # Sort the population (descending order)
        self._population.sort(reverse=True)

        # Remove _selection_rate of the population
        number_individuals = len(self._population)
        number_kept = round(self._selection_rate*number_individuals) # number of intividuals to keep, based on selection_rate
        # Remove less adapted individuals
        self._population = self._population[:number_kept]

        # Crossing (select 2 parents then create a new individual)
        for n in range(number_individuals-number_kept): # The quantity of new individuals is the same as the removed quantity
            # Apply the crossing function specific to the problem
            new_chromosome = self._problem.crossing_function(population=self._population, selection_rate=self._selection_rate) 
            # Generate a new individual from the above new_chromosome
            new_individual = Individual(chromosome=new_chromosome, fitness=self._problem.fitness_function(new_chromosome))
            # Add the new_individual to the population
            self._population.append(new_individual)

        # Apply the mutation function specific to the problem
        self._problem.mutation_function(self._population, self._mutation_rate)

        # Increment the generation count
        self.generation+=1 


    def showGenerationSummary(self):    
        """ Print some debug information on the current state of the population """
        # Display the current generation count
        print('Generation: ',self.generation) 
        # Display a header
        print('chromosome\tfitness')
        # Display the chromosome and fitness for each individuals of the population
        for i in self._population:
            print(i.chromosome,'\t',i.fitness)

    def getBestIndividual(self):
        """ Return the best Individual of the population """
        # Sort the population according to the fitness value (the higher the better)
        self._population.sort(reverse=True)
        # Get the first element of the population, which is the best individual
        best=self._population[0]
        return best

    def evolveUntil(self, max_nb_of_generations=500):
        """ Launch the evolveForOneGeneration function until one of the two condition is achieved : 
            - Max nb of generation is achieved
            - The fitness of the best Individual is greater than or equal to threshold_fitness
        """
        fitness_max = self.getBestIndividual().fitness # Maximal fitness of the population
        while (self.generation < max_nb_of_generations and fitness_max < self._problem._threshold_fitness):
            fitness_max = self.getBestIndividual().fitness
            self.evolveForOneGeneration()