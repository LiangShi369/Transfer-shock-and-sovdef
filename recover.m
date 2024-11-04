function aa = recover(a,x)

% construct full length sequence
aa = zeros(length(x),1);
dif = length(x) - length(a);

if dif >0

for i = 1:length(x)
    if x(i) == 1
        aa(i) = nan;
    elseif x(i) == 0 
        aa(i) = a(i-dif);
    end
end

else
    aa = a;
    
end

