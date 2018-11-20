function output = computeAH(As, Type)

    k = length(As);
    
    for i = 1:k
        invAs{i} = pinv(As{i});
    end
    
    Ar = As{1};
    for i = 2:k
       Ar = Ar + As{i};
    end
    Ar = Ar/k;
    
    if (Type == 0)
        output = Ar;
    end
    
    
    if(Type == 1)
    Ha = pinv(As{1});
    for i = 2:k
        Ha = Ha + invAs{i};
    end
    
    Ha = pinv(Ha/k);
    Arh = Ar^(0.5);
    Arinvh = Ar^(-0.5);
    temp = Arinvh * Ha * Arinvh;
    output = Arh * temp^(0.5) * Arh;
    end
end