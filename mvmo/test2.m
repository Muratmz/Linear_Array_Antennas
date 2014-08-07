%% MATLAB CODE FOR REAL CODED GENETIC ALGORITHM
%% Prepared by SOUMAVA KUMAR ROY,IIS2013001









clc;
clear all;
close all; 
display('Program to find the solution of f(x)=x^2 using genetic algorithm');
p_no=input('Enter the the size of initial population::');
bits=input('Enter the bitsize of each strings::') ;
final=[];
lower=input('Enter the lower limit of the variable x::');
upper=input('Enter the upper limit of the variable x::');
iterations=input('Enter the number of iterations to perform the process of genetic algorithm::');
tournament_no=input('Enter the number of times you want to perform tournament selection::');
p=rand(p_no,bits);
%% generate the bits randomly
for i=1:p_no
    for j=1:bits
        if (p(i,j)>=0.5)
            popu(i,j)=1;
        else 
            popu(i,j)=0;
        end;
        
    end;
end;
popu;
for z=1:iterations
val=zeros(p_no,1);
%%calculate the value of strings
for i=1:p_no
      for j=bits:-1:1
          k=bits-j;
          val(i,1)=val(i,1)+pow2(popu(i,j),k);
      end;
end;
val;
%%calculate the values mapped in the range [lower upper]
for i=1:p_no
   val_mapped(i,1)=lower+(((upper-lower)*val(i,1))/(pow2(bits)-1)); 
end
val_mapped;
%% generate the random numbers to do the tournament selection
for i=1:tournament_no
    select(i,:)=[floor(1+rand(1)*(p_no-1)) ceil(1+rand(1)*(p_no-1))];%% have to make changes here
    
end;
select;
%% calculate the fitness value and perform tournament selection
for i=1:tournament_no
    tour(i,:)=[val_mapped(select(i,1),1) val_mapped(select(i,2),1)];
    tour_val(i,:)=power(tour(i,:),2);
            if (tour_val(i,1)>=tour_val(i,2))
            fitest(i,:)=popu(select(i,1),:);
            else
            fitest(i,:)=popu(select(i,2),:);
            end;
end;
tour;
tour_val;
fitest;
prob_crossover=0.8;
prob_mutation=0.01;
crossover_point=2;
for x=1:ceil(p_no/2)
%%selecting of parents
parents=[];
crossover=[];
crossover_val=[];
crossover_val_mapped=[];
crossover_fitness=[];
parents=[fitest(floor(1+rand(1)*(tournament_no-1)),:); fitest(ceil(1+rand(1)*(tournament_no-1)),:)];
%%selecting of childrens
val_gen=rand(1);
child=[];
if (val_gen<prob_crossover)
   child=[[parents(1,1:crossover_point) parents(2,(crossover_point+1):bits)];[parents(2,1:crossover_point) parents(1,(crossover_point+1):bits)]];
end;
%%calculate the value of parents and children
crossover=[parents;child];
crossover_val=zeros(size(crossover,1),1);
for i=1:size(crossover,1)
      for j=bits:-1:1
          k=bits-j;
          crossover_val(i,1)=crossover_val(i,1)+pow2(crossover(i,j),k);
      end;
end;
crossover_val;
%%map the above value in the range [lower upper]
for i=1:size(crossover_val,1)
   crossover_val_mapped(i,1)=lower+(((upper-lower)*crossover_val(i,1))/(pow2(bits)-1)); 
end
crossover_val_mapped;
%% calculate the fitness value to select the the most suitable solution
crossover_fitness=power(crossover_val_mapped,2);
select_max1=max(crossover_fitness);
for i=1:size(crossover,1)
        if (select_max1==crossover_fitness(i,1))
            pos1=i;
            break;
        end;
end;
%%the position of the max fitness value after crossover
pos1;
%% to select the 2nd bst fitness value.
select_max2=crossover_fitness(1,:);
pos2=1;
if select_max2==select_max1
    select_max2=crossover_fitness(2,:);
    pos2=2;
end;
for i=1:size(crossover_fitness,1)
        if(select_max1==crossover_fitness(i,:))
            continue;
        elseif (crossover_fitness(i,:)>select_max2)
            select_max2=crossover_fitness(i,:);
            pos2=i;
        end;
end;
%% the values of 2nd best fitness values and its position
select_max2;
pos2;
%% taking the best 2 fitest solution and creating a different mating pool            
final=[final;crossover(pos1,:);crossover(pos2,:)];
end;%% end of crossover

%%performing mutation over the solutions generated after crossover.   
for i=1:size(final,1)
    for j=1:size(final,2)
    val_gen=rand(1);
        if (val_gen<prob_mutation)
          final1(i,j)=xor(final(i,j),1);
        else 
          final1(i,j)=final(i,j);
        end;
    end;
end;
final=[];  
popu=final1;
p_no=size(popu,1);
end;%%end to the number of iterations for the entire process;    
display('The final best solution is given as....')
popu