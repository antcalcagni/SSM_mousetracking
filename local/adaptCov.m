function [meann,C,k] = adaptCov(tempSample,meann,C,k,type_adapt)


switch type_adapt
    case 'type1'
        M=size(tempSample,2);
        xk = mean(tempSample,2);
        aux = tempSample - xk*ones(1,M);
        Ck = cov(aux',1);
        w1 = k/(k+1);
        w2 = 1/(k+1);
        w3 = k/(k+1)^2;
        meann = w1*meann + w2*xk;
        C = w1*C + w2*Ck + w3*(meann-xk)*(meann-xk)';
        C = topdm(C);
        k = k+1;
    case 'type2'
        xk = mean(tempSample,2);
        w2 = 1/(k+1);
        meann = meann + w2*(xk-meann);
        C = C + w2*( (xk-meann)*(xk-meann)' - C );
        C = topdm(C);
        k = k+1;
    case 'type3'
        %C = topdm(cov(tempSample'*2.4^2)/size(tempSample,1));
        C = topdm(cov(tempSample'));
        
end



end