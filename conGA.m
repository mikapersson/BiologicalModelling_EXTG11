function conGA
%exercise written 14 August 2006 by Anders Brodin

%global variables:
global POPSIZE GENE_NUMBER MAX_GENERATIONS NUM_BREED;
POPSIZE = 50;                   %popsize is number of chromosomes in the population
GENE_NUMBER = 2;                %gene number is number of genes in each chromosome, should be one in these simple tasks
MAX_GENERATIONS = 25;           %number of iterations
NUM_BREED = 20;                 %out of the 50 chromsomes (= possible solutions) the 20 best are allowed to reproduce,
                                %they double?

%genetic operators
global mutation_rate mutation_free;%crossover_rate 
mutation_rate = 0.04;           %how large percentage of the genes experience mutations
mutation_free = 10;             %number of topranked Individuals that should be saved from mutation

%global structure delaration
global Ind Gene Fitness;
Make_Array();
global Results1 Results2 Results3;
Results1 = zeros(MAX_GENERATIONS);
Results2 = zeros(MAX_GENERATIONS);
Results3 = zeros(MAX_GENERATIONS);

rand('twister',sum(100*clock));                     %'state'
Create_pop                                          %calls function that generates random chromosomes
for generation = 1 : MAX_GENERATIONS                %how long should the algorithm run?
   Evaluate_fitness;                                %calling function "evaluate_fitness" below
   Sort;  %highest fitness first (on top of Ind)    %calling function "sort" below
   LinReproduce;                                       %calling function "reproduce" below
   Mutate;                                          %calling function "mutate" below
   Output (generation);                             %calling function that writes result to Excel sheet
end                                                 %end of "for" loop
%Ind(1).Gene()
sol = Ind(1).Fitness;

%CALCULATE AVERAGE OF #RUNS
runs = 10;
res = [];
for i = 1 : runs
    Create_pop                                         
    for generation = 1 : MAX_GENERATIONS              
       Evaluate_fitness;                               
       Sort;  %highest fitness first (on top of Ind)   
       LinReproduce;                                      
       Mutate;                                          
       Output (generation);                             
    end 
    res = [res, Ind(1).Fitness];
end
res
sum(res)/runs

%__________________________________________________________________________
function Create_pop()                               %This function creates a random population of 50 chromosomes
global POPSIZE GENE_NUMBER Ind;                     % Gene Fitness;
for Count1 = 1 : POPSIZE                            %For counter for the chromosomes
    for Count2 = 1 : GENE_NUMBER                    %For counter for genes inside chromosomes
        g = rand(1) * 10;                           %g is a random number between 0 and 10
        Ind(Count1).Gene(Count2) = g;
    end
    Ind(Count1).Fitness = -1;                       % This fills the ind.Fitness varaibale with -1 to start with
end
%__________________________________________________________________________
function Evaluate_fitness()                         %This function calculates the base two number of ones and zeros
global POPSIZE Ind;  % Gene Fitness;
for Count1 = 1 : POPSIZE
   x = Ind(Count1).Gene(1);
   y = Ind(Count1).Gene(2);
   z = x*sin(4*x) + 1.1*y*sin(2*y);
   Ind(Count1).Fitness = z;
end
%__________________________________________________________________________
function Reproduce()
global POPSIZE NUM_BREED GENE_NUMBER Ind; % crossover_rate;
countback = POPSIZE;
for Count1 = 1 : 2 : (POPSIZE / 2)
    if Count1 <= NUM_BREED       
        Gene_No = 8 * rand(1) + 1;
        for Count2 = 1 : GENE_NUMBER
            if Count2 <= Gene_No
                Ind(countback).Gene(Count2) = Ind(Count1).Gene(Count2);
                Ind(countback - 1).Gene(Count2) = Ind(Count1 + 1).Gene(Count2);
            else  %vad innebär detta steget?
                Ind(countback).Gene(Count2) = Ind(Count1 + 1).Gene(Count2);
                Ind(countback - 1).Gene(Count2) = Ind(Count1).Gene(Count2);
            end
        end
        countback = countback - 2;        
    end
end
%_________
function PropReproduce()
global POPSIZE NUM_BREED GENE_NUMBER Ind; % crossover_rate;
countback = POPSIZE;
for Count1 = 1 : 2 : (POPSIZE / 2)
    if Count1 <= NUM_BREED       
        Gene_No = 8 * rand(1) + 1;
        for Count2 = 1 : GENE_NUMBER
            if Count2 <= Gene_No
                Ind(countback).Gene(Count2) = Ind(Count1).Gene(Count2);
                Ind(countback - 1).Gene(Count2) = Ind(Count1 + 1).Gene(Count2);
            else  %vad innebär detta steget?
                Ind(countback).Gene(Count2) = Ind(Count1 + 1).Gene(Count2);
                Ind(countback - 1).Gene(Count2) = Ind(Count1).Gene(Count2);
            end
        end
        countback = countback - 2;        
    end
end
%______________________________________________________________________________________________________________
function LinReproduce()
global POPSIZE NUM_BREED GENE_NUMBER Ind; % crossover_rate;
countback = POPSIZE;
for Count1 = 1 : 2 : (POPSIZE / 2)
    if Count1 <= NUM_BREED       
        p1 = Ind(Count1).Gene();
        p2 = Ind(Count1+1).Gene();
        offsp1 = 0.5*p1 + 0.5*p2;
        offsp2 = 1.5*p1 - 0.5*p2;
        offsp3 = 1.5*p2 - 0.5*p1;
        
        Ind(countback).Gene = mod(offsp1, 10);
        Ind(countback-1).Gene = mod(offsp2, 10);
        Ind(countback-2).Gene = mod(offsp3, 10);
        
        countback = countback - 3;        
    end
end
%______________________________________________________________________________________________________________
function Mutate()
global GENE_NUMBER POPSIZE mutation_free mutation_rate Ind;
for Count1 = mutation_free : POPSIZE  %shouldn't it be mutation_free + 1?
    for Count2 = 1 : GENE_NUMBER
        Chance = rand(1);
        if Chance <= mutation_rate
            mut = rand(1) * 10;
            Ind(Count1).Gene(Count2) = mut;
        end
    end
end
%__________________________________________________________________________
function Sort()
global POPSIZE Ind;
for i = 1 : (POPSIZE - 1)
    m = i + 1;  %index of gene with highest fitness
    for j = (i + 1) : POPSIZE  %why not just put m inst. of (i+1) here?
        if Ind(j).Fitness < Ind(m).Fitness
            m = j;
        end
    end
    if Ind(m).Fitness < Ind(i).Fitness  %finns det ingen swap funktion för Ind?
        tmp = Ind(i);
        Ind(i) = Ind(m);
        Ind(m) = tmp;
    end
end
%__________________________________________________________________________
function Output(generation)
global POPSIZE Ind MAX_GENERATIONS Results1 Results2 Results3;
addfitness = 0;
add_top = 0;
for Count = 1 : POPSIZE
    addfitness = addfitness + Ind(Count).Fitness;
    if Count <= 10
       add_top = add_top + Ind(Count).Fitness;
    end 
end
meantopfitness = add_top / 10;
meanfitness = addfitness / POPSIZE;
Results1(generation) = generation;
Results2(generation) = meanfitness;
Results3(generation) = meantopfitness;
%meantopfitness  %prints meantopfitness
if generation == MAX_GENERATIONS
    plot(Results1,Results2,Results1,Results3)
end
%________________________________________________________________________
function Make_Array()
global Ind Gene Fitness POPSIZE GENE_NUMBER;
for Count = 1 : POPSIZE
    Ind(Count).Gene = [];               
    Ind(Count).Fitness = -2;
    for Count2 = 1 : GENE_NUMBER
        Ind(Count).Gene(Count2)= 0;
    end
end
%____________________________________________________________





