function[c] = polyadd(p,q)
if length(p) > length(q)
    q=[zeros(1,length(p)-length(q)),q];
elseif length(q) > length(p)
    p=[zeros(q,length(q)-length(p)),q];
end

a=wrev(p); b=wrev(q);
y=a+b;
c=wrev(y);

    
