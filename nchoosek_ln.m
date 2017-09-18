function value = nchoosek_ln(a,b)
%This function takes the logarithm of nchoosek; needed to reduce numerical overflow. 
value = exp(log(factorial(a))-log(factorial(a-b))-log(factorial(b)));
end 