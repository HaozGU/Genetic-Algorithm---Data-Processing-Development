# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 2022

@author: tdrumond & agademer

Template file for your Exercise 3 submission 
(GA solving Mastermind example)
"""
from GASolver_module import GAProblem
import random
import mastermind_module as mm


class MastermindProblem(GAProblem):
    """Implementation of GAProblem for the mastermind problem"""

    def __init__(self, possible_genes, threshold_fitness):
        super().__init__(possible_genes, threshold_fitness)

    def generate_chromosome(self):
        chromosome = []
        # Each chromosome has size genes
        for g in range(size):
            chromosome.append(self._possible_genes[random.randint(0,len(self._possible_genes)-1)])
        return chromosome

    def fitness_function(self, chromosome):
        return match.rateGuess(chromosome)
        
    def crossing_function(self, population, selection_rate): # return new chromosome
        # Crossing the survivors to recreate the same quantity of individuals as originally
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
        new_chromosome = parent1[:2] + parent2[2:] # First two genes from parent1, last two genes from parent2
            
        return new_chromosome

    def mutation_function(self, population, mutation_rate):
        # Mutation for _mutation_rate of the population
        for i in population:
            number = random.random()
            if (number < mutation_rate):
                # A random value is assigned to a random gene of the chromosome
                i.chromosome[random.randint(0,3)] = self._possible_genes[random.randint(0,len(self._possible_genes)-1)]
                i.fitness = match.rateGuess(i.chromosome) # The fitness is recomputed for the mutated chromosome
        return 0


if __name__ == '__main__':

    from GASolver_module import GASolver

    size=6
    match = mm.MastermindMatch(secretSize=size)
    problem = MastermindProblem(possible_genes=list(range(len(mm.getPossibleColors()))), threshold_fitness=3*size)
    solver = GASolver(problem)

    solver.resetPopulation()
    solver.showGenerationSummary()
    print()
    solver.evolveUntil(5000)
    solver.showGenerationSummary()

    best = solver.getBestIndividual()

    print(best)
    print(f"Problem solved? {match.isCorrect(best.chromosome)}")

    #print(f"Best guess {mm.decodeGuess(solver.getBestIndividual().chromosome)} {solver.getBestIndividual()}")
    #print(f"Problem solved? {match.isCorrect(solver.getBestIndividual().chromosome)}")
