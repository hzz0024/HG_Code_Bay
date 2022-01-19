#### fitness(NULL) callback (manual p221): 
callback is not specific to a particular mutation type like a normal fitness() callback (see chapter 10), and it is not called on a per-mutation basis. Instead, fitness(NULL) callbacks are called once per individual per generation, and can provide a per-individual fitness effect that gets multiplied together with all of the other fitness effects in the model to produce the fitness of an individual (see section 25.2). In section 13.1, we will use a fitness(NULL) callback to implement the fitness function for a simple polygenic trait. 

#### fitnessScaling (manual p221): 
property is actually quite similar, but it is a property of Individual rather than being a callback. This property has a default value of 1.0, and it is reset to 1.0 in every generation, but script can set it to whatever fitness effect is desired for an individual. It value will be multiplied together with all of the other fitness effects in the model, just like the result of a fitness(NULL) callback, to produce the fitness of an individual. Using fitnessScaling can be faster and more convenient than a fitness(NULL) callback because it can be used in a vectorized fashion, allowing fitness effects to be determined across whole subpopulations with a single vectorized assignment. We will use fitnessScaling in most of the recipes in this chapter. 

#### make m2 colored:
m2.color = "red";

#### nonWF model non-overlapping generations; kill off the parental generation (p402)     
inds = sim.subpopulations.individuals;
// non-overlapping generations; kill off the parental generation
ages = inds.age;
inds[ages > 0].fitnessScaling = 0.0;
inds = inds[ages == 0];

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