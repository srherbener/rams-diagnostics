function fact = myFact( n )
%myFact return factorial of n
%   factorial 0 --> 1
%   factorial 1 --> 1

fact = 1;
for i = 2:n
    fact = fact * i;
end


end

