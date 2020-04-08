function TSP

%global variables:
global POPSIZE GENE_NUMBER MAX_GENERATIONS NUM_BREED DISTANCES;
POPSIZE = 50;                   %popsize is number of chromosomes in the population
GENE_NUMBER = 5;                %gene number is number of genes in each chromosome, should be one in these simple tasks
MAX_GENERATIONS = 10;           %number of iterations
NUM_BREED = 20;                 %out of the 50 chromsomes (= possible solutions) the 20 best are allowed to reproduce
DISTANCES = [0,267,602,506,262,282;     %distances as given in exercise Q10
             267,0,412,467,341,214;
             602,412,0,305,467,321;
             506,467,305,0,247,243;
             262,341,467,247,0,149;
             282,214,321,243,149,0];

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

rand('twister',sum(100*clock));%'state'
Create_pop                                          %calls function that generates random chromosomes
for generation = 1 : MAX_GENERATIONS                %how long should the algorithm run?
   Evaluate_fitness;                                %calling function "evaluate_fitness" below
   Sort;  %highest fitness first (on top of Ind)    %calling function "sort" below
   Reproduce;                                       %calling function "reproduce" below
   Mutate;                                          %calling function "mutate" below
   Output (generation);                             %calling function that writes result to Excel sheet
end                                                 %end of "for" loop
Ind(1).Gene
%__________________________________________________________________________
function Create_pop()                               %This function creates a random population of 50 chromosomes
global POPSIZE GENE_NUMBER Ind;                     % Gene Fitness;
for Count1 = 1 : POPSIZE                            %For counter for the chromosomes
    for Count2 = 1 : GENE_NUMBER                    %For counter for genes inside chromosomes
        g = floor(2 + GENE_NUMBER*rand(1));                          %g is a random number between 0 and 1
        if g == 7
            g = 6;
        end
        
        exists = sum(ismember(Ind(Count1).Gene(), g));  %check if city exists in route/individual
        while exists ~= 0
            g = floor(2 + 5*rand(1));                          
            if g == 7
                g = 6;
            end
            exists = sum(ismember(Ind(Count1).Gene(), g));
        end
        
        Ind(Count1).Gene(Count2) = g;
    end
    Ind(Count1).Fitness = -1;                       % This fills the ind.Fitness varaibale with -1 to start with
end
%__________________________________________________________________________
function Evaluate_fitness()                         %This function calculates the route distance
global POPSIZE GENE_NUMBER Ind DISTANCES;% Gene Fitness;
for Count1 = 1 : POPSIZE
   Added_Value = DISTANCES(1, Ind(Count1).Gene(1));
   for Count2 = 1 : GENE_NUMBER-1
       dist = DISTANCES(Ind(Count1).Gene(Count2),Ind(Count1).Gene(Count2+1));
       Added_Value = Added_Value + dist;
   end
   Added_Value = Added_Value + DISTANCES(1,Ind(Count1).Gene(5));
   Ind(Count1).Fitness = Added_Value;
end
%__________________________________________________________________________
function Reproduce()
global POPSIZE NUM_BREED GENE_NUMBER Ind;   %two parents create one child through ordered crossover
countback = POPSIZE;
for Count1 = 1 : 2 : (POPSIZE / 2)
    if Count1 <= NUM_BREED       
        child = zeros(1,GENE_NUMBER);
        geneA = floor(rand(1)*GENE_NUMBER + 1);
        geneB = floor(rand(1)*GENE_NUMBER + 1);
        
        if geneA == 6
            geneA = 5;
        end
        if geneB == 6
            geneB = 5;
        end
        
        startGene = min(geneA,geneB);
        endGene = max(geneA,geneB);
        
        for i = startGene : endGene
            child(i) = Ind(Count1).Gene(i);
        end
        
        Count2 = 1;
        while Count2 < 6
           if child(Count2) == 0
              Count3 = 1;
              exists = sum(ismember(child, Ind(Count1+1).Gene(Count3)));  %avoid multiples in child
              while exists ~= 0
                  Count3 = Count3+1;
                  exists = sum(ismember(child, Ind(Count1+1).Gene(Count3)));
              end               
              child(Count2) = Ind(Count1+1).Gene(Count3);
           end
           Count2 = Count2 + 1;
        end
        Ind(countback).Gene = child;
        countback = countback - 1;        
    end
end
%__________________________________________________________________________
function Mutate()
global GENE_NUMBER POPSIZE mutation_free mutation_rate Ind;
for Count1 = mutation_free+1 : POPSIZE
    Chance = rand(1);
    if Chance <= mutation_rate
        g = floor(rand(1)*GENE_NUMBER+1);
        if g == 6
            g = 5;
        end
        Ind(Count1).Gene([g, (GENE_NUMBER+1-g)]) = Ind(Count1).Gene([(GENE_NUMBER+1-g), g]);
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
    for Count2 = 1 : GENE_NUMBER
        Ind(Count).Gene(Count2)= 0;
    end
end





