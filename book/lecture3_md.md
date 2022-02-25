# Mutation, Migration, with and without selection

# Mutation and Evolution 

Previously we learned about how selection can change the frequencies of
alleles and genotypes in populations. Selection typically eliminates
variation from within populations. (The general exception to this claim
is with the class selection models we have called "balancing\" selection
where alleles are maintained in the population by overdominance,
habitat-specific selection, or frequency dependent selection). If
selection removes variation, soon there will be no more variation for
selection to act on, and evolution will grind to a halt, right? This
might be true if it were not for the reality of mutation which will
restore genetic variation eliminated by selection. Thus, mutations are
the fundamental raw material of evolution.

The basic model of mutation that we will study is one way mutation. This
is when one allele through mutation can turn in to another such that

$$\begin{aligned}
A_1 \stackrel{u}\longrightarrow A_2
\end{aligned}$$

$u$ here represents the mutation rate, the probability that a mutation
from $A_1$ to $A_2$ occurs during a meiosis.

It should be obvious that mutation will change allele frequencies. This
is true because there is a constant flux from $A_1$ to $A_2$ purely as a
result of this mutation process. We can study the change in allele
frequency due to mutation in a very similar way to how we studied the
change in allele frequency due to selection. Consider a population with
frequency of the $A_1$ allele $p$. In the next generation, after a round
of mutation, each $A_1$ allele must have been $A_1$ in the current
generation and it must not have mutated. That is 

$$\begin{aligned}
    p' = p (1-u).
\end{aligned}$$ 

Now lets turn our attention to the
change in allele frequency in one generation due to mutation as we did
previously for selection 

$$\begin{aligned}
    \Delta_up & = & p' - p \\
     & = & p(1-u) -p \\
     & = & -up \\
\end{aligned}$$

 notice that our notation-- $\Delta_up$
-- emphasizes the source of the change in allele frequency is mutation.

This sort of unidirectional mutation acts to consistently decrease the
frequency of the $A_1$ allele from generation to generation. If we
instead were to study two way mutation, the mutational flux would depend
on the proportional rates to and from $A_1$. So mutation, although it is
a random process with respect to target, leads to deterministic effects
on allele frequencies. Neat huh?

Mutation rates per generation are very very small. You know this
intuitively- think about cloning plants from cuttings. In *Drosophila*,
which is one of the best studied animals from the perspective of rates
of spontaneous mutation, the mutation rate per generation per nucleotide
is on the order of $10^{-9}$. This means that mutation changes the
frequency of alleles are a very slow rate. If we assume that there is
sufficiently strong selection against the $A_2$ allele (i.e.
$w_{11} >> w_{22}$), then we can further approximate our change in
allele frequency due to mutation 

$$\begin{aligned}
\Delta_up  & = & -up \\
  & = & -u + qu \\
 & \approx & -u, 
\end{aligned}$$ 

because $q \approx 0$. So in the case
of a deleterious $A_2$ allele, we can see that the change in allele
frequency due to mutation is independent of allele frequencies. This is
our first hint that mutation and selection might combine in interesting
and important ways.

# Mutation-Selection Balance 

We can imagine mutation and selection as opposing forces which might
come to some equilibrium in terms of the number or frequency of
deleterious alleles within a population. To make this concrete think
about the human genetic disease cystic fibrosis (CF). CF is a very
serious genetic disorder in which a transmembrane protein in lung
epithelium cells called CFTR is non-functional. Hundreds, if not
thousands of separate mutations in CFTR lead to CF, thus we could
imagine that there is a certain, appreciable rate of mutation to CF. If
each of these mutations is deleterious (i.e. they cause disease) then
over generations they should be selected out of the population. Thus
mutation will inject CF mutations into the population, but selection
will remove them- can we study this as an equilibrium process?

Our approach will be to study each of our evolutionary forces in
isolation, and then combine them to figure out how they interact. Let's
start by considering the change in allele frequency due to selection
that we studied in lecture 7, but this time we will approximate it under
the assumption that $q \approx 0$ 

$$\begin{aligned}
\Delta_sp  & = &\frac{pqs[ph + q(1-h)]}{\bar{w}} \\
& \approx & qhs 
\end{aligned}$$ 

This approximation goes down the road
because when $q \approx 0$, $p \approx 1$, $\bar{w} \approx 1$, and we
can ignore all terms of order $q^2$.

Now lets combine the forces of selection and mutation on the change in
frequency of $A_1$ using the approximations we have just derived
(equations 3 and 5). At equilibrium that change in allele frequency due
to the combined actions of mutation and selection must equal zero. That
is 

$$\begin{aligned}
    0 & = & \Delta_up + \Delta_sp \\
    & \approx & -u + qhs
\end{aligned}$$ 

so the equilibrium frequency of
the $A_2$ is 

$$\begin{aligned}
    \hat{q} \approx \frac{u}{hs}
\end{aligned}$$

Thus we see that
deleterious (e.g. disease) allele frequencies are determined by both the
mutation rate to those alleles and their selective effects in
heterozygotes. As we saw earlier, new mutations overwhelmingly are found
in heterozygous states, so it's perhaps not surprising that $h$ should
dominate the fate of deleterious alleles.

**Problem:** Go through the same steps of approximations that we just
did to find the Mutation-Selection equilibrium value of mutations which
are completely recessive (i.e. $h = 0$). This would make a heck of an
exam question\....

# Migration and Gene Flow 

In population genetics, the term "migration\" is really meant to
describe Gene flow, defined as the movement of alleles from one area
(deme, population, region) to another. Gene flow assumes some form of
dispersal or migration (wind pollination, seed dispersal, birds flying,
etc.) but dispersal is not gene flow (genes must be transferred, not
just their carriers)

We are going to build a model of gene flow in exactly the same way we
studied mutation. Consider two populations, a mainland population and an
island population. Each of these populations has the $A_1$ allele at
frequencies $p_{main}$ and $p_{island}$ respectively. Assume that gene
flow is one way, from mainland to island and that the proportion of
individuals who become parents in the island population is $m$. Although
I've said this is a proportion, also notice that we could consider this
a probability interchangeably. After a round of migration, in the island
population there are then two sources for alleles, they could be from
the island population originally with probability $1-m$, or they could
have migrated from the mainland with probability $m$. This means that
after migration the allele frequency of $A_1$ in the island population
is 
$$\begin{aligned}
p_{island}' = p_{island}(1-m) + p_{main}m.
\end{aligned}$$ 

Simple enough right? Now lets to what comes naturally and study the change in allele
frequency as a result of mutation. Following what we have done in our
other analyses, 

$$\begin{aligned}
\Delta_{m}p_{island}  & = & p_{island}' - p_{island} \\
 & = & p_{island}(1-m) + p_{main}m - p_{island} \\
 & = & m(p_{main} - p_{island}) \\
\end{aligned}$$ 

Beautiful. Now we have
a very simple expression for how allele frequencies in the island
population should change due to gene flow from the mainland population.
This change in allele frequency makes sense- it only depends on the
amount of migration and the differences in allele frequencies between
the two populations. It's simple enough to generalize this to multiple
populations, or to populations at different distances away from one
another, but we won't cover that here.

