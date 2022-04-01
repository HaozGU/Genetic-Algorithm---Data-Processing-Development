# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 2022

@author: tdrumond

Module containing utility functions to implement a Mastermind game.

The main element is the MastermindMatch class, that allows to instantiate 
a match of the game having a certain secret code generated at random.
This class plays the role of the codemaker, allowing to check if a guess
is correct and rating how close a guess is to the secret code.
"""
from random import choice
from typing import List

# Possible colors for codes in in the game
_colors = ['blue', 'red', 'green', 'yellow', 'orange', 'violet']


def getPossibleColors():
    """Getter function to read the array of possible colors"""
    return _colors


def decodeGuess(guess: List[int]) -> List[str]:
    """Decode a guess from integers to their corresponding colors

    Args:
        guess (list[int]): a mastermind guess as a list of integers

    Returns:
        list[str]: a mastermind guess as a list of color strings
    """
    return [_colors[i] for i in guess]


def generateRandomSecret(size):
    """Generate a random secret of a given size"""
    secret = [choice(range(len(_colors))) for i in range(size)]
    return secret


class MastermindMatch:
    """Class to instantiate a mastermind game match with a random secret code.
    A MastermindMatch object plays the role of the codemaker player, while the class user typically plays the role of the code guesser.
    """

    def __init__(self, secretSize=4, correctColorPoints=1, correctPositionPoints=3):
        """Instantiates a mastermind guess with a random secret code

        A match can be created by calling:
        match = MastermindMatch()

        Args:
            secretSize (int, optional): defines the size of the secred. Defaults to 4.
            correctColorPoints (int, optional): points awarded for a correct color at the wrong position. Defaults to 1.
            correctPositionPoints (int, optional): points awarded for a correct color at the right position. Defaults to 3.
        """
        self._secret = generateRandomSecret(secretSize)
        self.correctColorPoints = correctColorPoints
        self.correctPositionPoints = correctPositionPoints

    def isCorrect(self, guess: List[int]) -> bool:
        """Checks whether a guess matches the secret code

        Args:
            guess (list[int]): a mastermind guess as a list of integers

        Returns:
            bool: True if the guess matches the secret code, False otherwise
        """
        return guess == self._secret

    def generateRandomGuess(self):
        return generateRandomSecret(len(self._secret))

    def rateGuess(self, guess: List[int]):
        """Gives a numeric score for a given guess proportional to how close 
        it is to the secret code (higher is better)

        Args:
            guess (list[int]): a mastermind guess as a list of integers

        Returns:
            int or float: the computed score
        """
        
        correctPosition = 0
        correctColors = 0
        print(guess)
        for i, color in enumerate(guess):
            if self._secret[i] == color:
                correctPosition += 1
            elif color in self._secret:
                correctColors += 1
        score = correctColors*self.correctColorPoints + \
            correctPosition * self.correctPositionPoints
        #print("Comparing ", guess, "with", self._secret, ":",score)
        return score

    def secretSize(self):
        """Returns the size of the secret code"""
        return len(self._secret)

    def maxScore(self):
        """Returns the maximum possible score under the defined point 
        schedule for this instance"""
        return self.correctPositionPoints * len(self._secret)
