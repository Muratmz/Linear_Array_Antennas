%Matyas function
%supports vectorized input
%pop is the 1 X nVar sized matrix
%in case of vectorized input, it is vctr_size X nVar matrix
function val = matyas(pop)
val = zeros(size(pop,1),1);
for ii = 1:size(pop,1)
    p = pop(ii,:);
    val(ii,:) = 0.26*(p(1).^2 + p(2).^2)-0.48*p(1).*p(2);
end

