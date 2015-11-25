function [ new_a ] = poole_stable(a)

[R,P,K] = residue(1,a);
for i=1:length(P)
    dist = norm(P(i));
    if dist > 1
        P(i) = 1/P(i);
    end  
end

[b,new_a] = residue(R,P,K);

end

