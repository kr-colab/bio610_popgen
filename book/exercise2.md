---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---


# Exercise 2

For our in class exercise you programmed a little
function to compute the allele frequency in the next 
generation, $p'$, due to natural selection.
I'd like you to use this function to explore the 
_sojourn time_ of a selected mutation in the population.
The sojourn time is the amount of time it takes for the allele
to go from some initial frequency, $p_0$, to frequency 1 (fixation) in a population. 

For this exercise, make a plot of sojourn time as a function of
the selection coefficient $s$ for perhaps 100 values of $s$
from between zero and one, with a starting frequency $p_0 = 0.001$.