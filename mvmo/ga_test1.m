%find the global minima of the rastrigen function using genetic algoithm
Fitnessfunction = @rastriginsfcn;
Noofvar = 2;
[x,fval] = ga(Fitnessfunction,Noofvar)

disp 'testing';
disp '';
rastriginsfcn(0) 
disp '';
size(x) 