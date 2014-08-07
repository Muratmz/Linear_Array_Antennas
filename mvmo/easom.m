%easom function
function val = easom(pop)
val = -1*cos(pop(1))*cos(pop(2))*exp(-1*((pop(1) - pi)^2 + (pop(2) - pi)^2));