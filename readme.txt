- NSGA2.m : The main code of NSGA-II
  - Calculate_Pareto_Front.m : Obtain the ideal pareto front
  - Pareto_Set.m : Get all individuals in a front
  - Pareto_Rank.m: Calculate the pareto rank of individuals
  - Crowding_Distance: Calculate the crowding distance of individuals
- MOEA_D.m : The main code of MOEA/D
  - UniformPoint.m : Obtain a set of uniform weight vectors
- drawFig.m : Plot the convergencing process

DTLZ test suite:
  - DTLZ1.m
  - DTLZ2.m
  - DTLZ3.m
  - DTLZ4.m
  - DTLZ5.m
  - DTLZ6.m
  - DTLZ7.m

- note.m : Execute the main algorithm and collect the result. Write all results in a excel file and in three different sheets

- optiFig.m : Plot the pareto front with minimum IGD value

- significance.m : Read the file and calculate the corresponding significance level

Besides, all the recorded results that shown in this project have been collected in File named 'IGD'
- DTLZx.xlsx
- NSGA2 : all the collected result in seven test problems obtained by NSGA-II with 2-objective and 3-objective 
- MOEA/D : all the collected result in seven test problems obtained by MOEA/D with 2-objective and 3-objective 