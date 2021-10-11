#### fitness(NULL) callback (manual p221): 
callback is not specific to a particular mutation type like a normal fitness() callback (see chapter 10), and it is not called on a per-mutation basis. Instead, fitness(NULL) callbacks are called once per individual per generation, and can provide a per-individual fitness effect that gets multiplied together with all of the other fitness effects in the model to produce the fitness of an individual (see section 25.2). In section 13.1, we will use a fitness(NULL) callback to implement the fitness function for a simple polygenic trait. 

#### fitnessScaling (manual p221): 
property is actually quite similar, but it is a property of Individual rather than being a callback. This property has a default value of 1.0, and it is reset to 1.0 in every generation, but script can set it to whatever fitness effect is desired for an individual. It value will be multiplied together with all of the other fitness effects in the model, just like the result of a fitness(NULL) callback, to produce the fitness of an individual. Using fitnessScaling can be faster and more convenient than a fitness(NULL) callback because it can be used in a vectorized fashion, allowing fitness effects to be determined across whole subpopulations with a single vectorized assignment. We will use fitnessScaling in most of the recipes in this chapter. 

#### make m2 colored:
m2.color = "red";


 