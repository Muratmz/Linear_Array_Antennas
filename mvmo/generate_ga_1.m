%function [x,fval,exitflag,output,population,score] = generate_ga_1
%    profile on
tic;

global g_t_thisPopulation 
Fitnessfunc = @booth;
%Fitnessfunc = @easom;

%% Start with the default options
options = gaoptimset;
%% Modify options setting
%options = gaoptimset(options,'MutationFcn', @mutfun);
options = gaoptimset(options,'MutationFcn', @mutationgaussian);
options = gaoptimset(options,'Display', 'off');
%options = gaoptimset(options,'PopInitRange', [-0.1; 0.1]);
options = gaoptimset(options,'Generations',100);
%options = gaoptimset(options,'PopulationSize',10);
%options = gaoptimset(options,'EliteCount',2);

[x,fval,exitflag,output,population,score] = ...
ga(Fitnessfunc,2,[],[],[],[],[],[],[],[],options);
fval 
x
% profile off
% profview
%Fitnessfunc = @roenbrock;
%Fitnessfunc = @matyas;
%Fitnessfunc = @rastriginsfcn;

toc
