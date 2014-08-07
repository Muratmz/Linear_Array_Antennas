%function [x,fval,exitflag,output,population,score] = generate_ga_1
%    profile on

global g_t_thisPopulation 
Fitnessfunc = @mvmo_fitness;
%Fitnessfunc = @easom;

%% Start with the default options
options = gaoptimset;
%% Modify options setting
%options = gaoptimset(options,'MutationFcn', @mutfun);
options = gaoptimset(options,'MutationFcn', @mutationgaussian);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'TolFun',1e-15);
options = gaoptimset(options,'PopInitRange',[-pi ;pi]);
options = gaoptimset(options,'Generations',100);
options = gaoptimset(options,'PopulationSize',10);
%options = gaoptimset(options,'EliteCount',2);

[x,fval,exitflag,output,population,score] = ...
ga(Fitnessfunc,32,[],[],[],[],[],[],[],[],options);

fval

% profile off
% profview
%Fitnessfunc = @roenbrock;
%Fitnessfunc = @matyas;
%Fitnessfunc = @rastriginsfcn;

