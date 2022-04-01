# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 2022

@author: agademer

Module containing utility functions to instantiate a traveling 
salesperson problem.
This module is accompanied of a cities.txt file, containing a list of
2D coordinates representing different cities.
"""

import matplotlib.pyplot as plt
from random import shuffle

def countCities(filename):
    """ return the number of cities in the file"""
    with open(filename) as file:
        nbCities = int(file.readline())
        return nbCities


def loadCities(filename):
    """ load the cities list from a text file """
    with open(filename) as file:
        cities = []
        next(file)
        for _ in range(countCities(filename)):
            x, y = file.readline().split(" ")
            cities.append((int(x), int(y)))
        return cities


def defaultRoad(cities):
    """ Default road: all the cities in the order of the text file """
    return list(range(len(cities)))


def drawCities(cities, road=None):
    """ Plot the cities and the trajectory """
    citiesCoordinates = list(zip(*cities))
    plt.figure()
    plt.scatter(citiesCoordinates[0], citiesCoordinates[1], color="red")
    if road is not None:
        roadCities = [cities[c] for c in road]
        # print(roadCities)
        roadCitiesCoordinates = list(zip(*roadCities))
        plt.plot(roadCitiesCoordinates[0], roadCitiesCoordinates[1])
        for i in range(len(roadCities)):
            plt.annotate(i, roadCities[i], xytext=(
                3, 3), textcoords='offset points')
        plt.gca().set_aspect('equal')
    plt.show


def distance(city1, city2):
    """ Euclidian distance between two cities """
    return ((city1[0] - city2[0])**2 + (city1[1] - city2[1])**2)**0.5


def roadLength(cities, road):
    """ Calculate the length of the road """
    roadCities = [cities[c] for c in road]
    total = 0
    for i in range(len(roadCities)-1):
        total += distance(roadCities[i], roadCities[i+1])
    total += distance(roadCities[-1], roadCities[0])
    return total


if __name__ == '__main__':
    cities = loadCities("cities.txt")
    print(cities)
    road = defaultRoad(cities)
    # road.reverse()
    shuffle(road)
    print(road)
    drawCities(cities, road)
    print(roadLength(cities, road))
