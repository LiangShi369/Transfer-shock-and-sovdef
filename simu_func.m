function [zpath] = simu_func(pdfz,T)

nz = size(pdfz,1) ;
zpath = zeros(T,1) ; 
zpath(1,1) = round(nz/2) ;  %start zpath at z closest to muz
    
for t = 2:T
    
    z0 = zpath(t-1,1) ;
    temp = cumsum( pdfz(z0,:) );
    shock = rand ;

    if shock < temp(1)
        zpath(t,1) = 1;
    elseif shock > temp(nz)
        zpath(t,1) = nz;
    else
        zpath(t,1) = find(shock<=temp(2:nz) & shock>temp(1:nz-1)) + 1; 
    end

end


end

