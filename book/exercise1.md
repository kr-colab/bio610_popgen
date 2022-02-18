# Exercise 1

Imagine that you and your friends have collected a number of [ice crawlers](https://en.wikipedia.org/wiki/Grylloblattidae) from Mt. Hood
and are interested in examining variation at the CLD1 locus, which is responsible for cold tolerance in this incredible insect. You find
two alleles at the CLD1 locus, the $C_1$ and the $C_2$ alleles, and in genotyping 50 individuals you find the following numbers
of genotypes.


| $C_1C_1$ | $C_1C_2$ | $C_2C_2$ |
|----:|----:|----:|
| 7   |  35  |  8 |


1.  Are these genotypes close to the proportions we would expect
    under Hardy-Weinberg? To answer this write a python function that takes as 
    input the three genotype counts above, computes the expected HW 
    proportions, and then performs a $\Chi^2$
    test using the `scipy` function
    [`scipy.stats.chisquare()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chisquare.html). 

    Ideally your function would return the expected HW proportions along with 
    the results of the statistical test. You could do this by returning a tuple...