#### fitness(NULL) callback (manual p221): 
callback is not specific to a particular mutation type like a normal fitness() callback (see chapter 10), and it is not called on a per-mutation basis. Instead, fitness(NULL) callbacks are called once per individual per generation, and can provide a per-individual fitness effect that gets multiplied together with all of the other fitness effects in the model to produce the fitness of an individual (see section 25.2). In section 13.1, we will use a fitness(NULL) callback to implement the fitness function for a simple polygenic trait. 

#### fitnessScaling (manual p221): 
property is actually quite similar, but it is a property of Individual rather than being a callback. This property has a default value of 1.0, and it is reset to 1.0 in every generation, but script can set it to whatever fitness effect is desired for an individual. It value will be multiplied together with all of the other fitness effects in the model, just like the result of a fitness(NULL) callback, to produce the fitness of an individual. Using fitnessScaling can be faster and more convenient than a fitness(NULL) callback because it can be used in a vectorized fashion, allowing fitness effects to be determined across whole subpopulations with a single vectorized assignment. We will use fitnessScaling in most of the recipes in this chapter. 

#### make m2 colored:
m2.color = "red";

#### nonWF model non-overlapping generations; kill off the parental generation (p402)     

```slim
inds = sim.subpopulations.individuals;
// non-overlapping generations; kill off the parental generation
ages = inds.age;
inds[ages > 0].fitnessScaling = 0.0;
inds = inds[ages == 0];
```

#### set up the migration rate of individuals to adjacent subpops in nonWF model (p402)

```Slim
numMigrants = rbinom(1, inds.size(), M);  
if (numMigrants)
{
migrants = sample(inds, numMigrants);
currentSubpopID = migrants.subpopulation.id;
displacement = -1 + rbinom(migrants.size(), 1, 0.5) * 2; // -1 or +1
newSubpopID = currentSubpopID + displacement;
actuallyMoving = (newSubpopID >= 0) & (newSubpopID < N);
if (sum(actuallyMoving))
{
migrants = migrants[actuallyMoving];
newSubpopID = newSubpopID[actuallyMoving];
// do the pre-planned moves into each subpop in bulk
for (subpop in sim.subpopulations)
subpop.takeMigrants(migrants[newSubpopID == subpop.id]);
}
}
```

#### sim.addSubpop can be written in two ways

One is

```Slim
1 { sim.addSubpop("p1", 500); }
```

Another one is

```Slim
1 early() {
	sim.addSubpop("p1", 500);
}
```

#### 16.7 Evolutionary rescue after environmental change

```slim
early() {
// QTL-based fitness
	inds = sim.subpopulations.individuals;
	phenotypes = inds.sumOfMutationsOfType(m1);
	optimum = (sim.generation < Tdelta) ? opt1 else opt2;
	deviations = optimum - phenotypes;
	fitnessFunctionMax = dnorm(0.0, 0.0, 5.0);
	adaptation = dnorm(deviations, 0.0, 5.0) / fitnessFunctionMax;
	inds.fitnessScaling = 0.1 + adaptation * 0.9;
	inds.tagF = phenotypes; // just for output below
	// density-dependence with a maximum benefit at low density
	p1.fitnessScaling = min(K / p1.individualCount, 1.5);
}
fitness(m1) { return 1.0; }
```

#### Interaction types

There are currently four options for interaction functions (IFs) in SLiM, represented by singlecharacter codes:

"f" – a fixed interaction strength. This IF type has a single parameter, the interaction strength to be used for all interactions of this type. By default, interaction types use a type "f" IF with a value of 1.0, so interactions are binary: on within the maximum distance, off outside.

"l" – a linear interaction strength. This IF type has a single parameter, the maximum interaction strength to be used at distance 0.0. The interaction strength falls off linearly, reaching exactly zero at the maximum distance. In other words, for distance d, maximum interaction distance dmax, and maximum interaction strength fmax, the formula for this IF is f(d) = fmax(1 − d / dmax). 

"e" – A negative exponential interaction strength. This IF type is specified by two parameters, a maximum interaction strength and a shape parameter. The interaction strength falls off nonlinearly from the maximum, and cuts off discontinuously at the maximum distance; typically a maximum distance is chosen such that the interaction strength at that distance is very small anyway. The IF for this type is f(d) = fmaxexp(−λd), where λ is the specified shape parameter. Note that this parameterization is not the same as for the Eidos function rexp().

"n" – A normal interaction strength (i.e., Gaussian, but "g" is avoided to prevent confusion with the gamma-function option provided for, e.g., MutationType). The interaction strength falls off 

"c" – A Cauchy-distributed interaction strength. The interaction strength falls off non-linearly from the maximum, and cuts off discontinuously at the maximum distance; typically a maximum distance is chosen such that the interaction strength at that distance is very small anyway. This IF type is specified by two parameters, a maximum interaction strength and a scale parameter. The IF for this type is f(d) = fmax/(1+(d/λ)2), where λ is the scale parameter. Note that this parameterization is not the same as for the Eidos function rcauchy(). A Cauchy distribution can be used to model interactions with relatively fat tails. 

#### The fitness(NULL) callback 

The fitness(NULL) callback sets up the conditions for polygenic selection. It is called once per
individual per generation, and implements the fitness function that converts individual phenotypic
trait values into fitness effects.

#### Difference between WF and non-WF model

1. In WF models, early() events occur just before offspring generation and late() events happen just after. In nonWF models this positioning is different, because offspring generation has been moved earlier in the generation cycle – now, first() events occur just before offspring generation and early() events happen just after.

2. in WF models late() events are usually the best place to add new mutations, because the next stage is fitness evaluation, which immediately incorporates the fitness effects of the new mutations; but in nonWF models early() events are usually the best place to add new mutations, for the same reason.

3. In WF models, late() events are generally the best place to put output events so that they reflect the final state at the end of a generation; in nonWF models, however, output can make sense in either an early() or a late() event, depending upon whether you want to see the state of the population before or after viability selection.

#### Viability selection in nonWF models 

Viability selection in nonWF models is mechanistically simple. For a given individual, a fitness
of 0.0 or less results in certain death; that individual is immediately removed from its
subpopulation. A fitness of 1.0 or greater results in certain survival; that individual remains in its
subpopulation, and will live into the next generation cycle. A fitness greater than 0.0 and less than
1.0 is interpreted as a survival probability; SLiM will do a random draw to determine whether the
individual survives or not.

#### Basic elements

1. In SLiMGUI, each yellow square represents one individual; there are 500 squares since the subpopulation size is 500. Each individual is colored according to the calculated fitness of that individual;

2. If all mutations are neutral, then all individuals have the same relative fitness, 1.0, and so they are all yellow (as expressed by the color stripe for fitness values, discussed in section 3.1). In a simulation with beneficial and deleterious mutations you would see a range of colors, giving you an immediate visual sense of the fitness distribution across the subpopulation

3. 1000: {... omitting the end generation. Just as supplying no generation range at all for an event means “run in every generation”, omitting the end generation means “run in every generation after the specified start generation”. One may similarly omit the start generation with the syntax :end, meaning “run in every generation from the beginning until the specified end generation”.

4. addSubpop() sets the age of all new individuals to 0, that would provide our model with a bit of an artificial start

5. The fitness(m1) callback makes the direct fitness effect of m1 mutations neutral, as usual in such QTL models; the only effect of QTL mutations on fitness is indirect, through their effect on individual phenotypic values (see chapter 13 for an introduction to QTL models).