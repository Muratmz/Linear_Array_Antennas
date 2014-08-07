function mc = submutfun(inp,parents,GenomeLength)
m = mean(inp);
v = var(inp);
mutationChildren = zeros(length(parents),GenomeLength);
fs = 0.2;
% 
%  for i=1:length(parents)
%         parent = thisPopulation(parents(i),:);
%         mutationChildren(i,:) = parent  + scale .* randn(1,length(parent));
%  end
s = -1*log(v)*fs;


s1 = 1;
s2 = 1;
ran= rand(size(mutationChildren,1),size(mutationChildren,2));
% 
% for cc = 1:size(mutationChildren,2)
%     
%     
% end



for rr = 1 :size(mutationChildren,1)
    for cc = 1:size(mutationChildren,2)
        
    h(rr,cc) =   m(cc)*(1-exp(-1*ran(rr,cc)*s1)) + (1-m(cc))*exp(1-ran(rr,cc)*s2);
    h0(cc) =     m(cc)*(1-exp(-1*0*s1)) + (1-m(cc))*exp(1-0*s2);
    h1(cc) =     m(cc)*(1-exp(-1*1*s1)) + (1-m(cc))*exp(1-1*s2);
        
        
    end
end

for rr = 1 :size(mutationChildren,1)
    for cc = 1:size(mutationChildren,2)

        mc(rr,cc) = h(rr,cc)+(1-h1(cc)+h0(cc))*ran(rr,cc)-h0(cc);
    end
end

%testing the part of code implemented in mutationgaussian.m by adding the
%mutated values got my my algorithm to the parents selected and passed to
%this function
for i=1:length(parents)
        parent = inp(parents(i),:);
        mc1(i,:) = parent  + mc(i,:);
        
end
mc = mc1;
%fprintf('mc value is :\n');
%mc;



