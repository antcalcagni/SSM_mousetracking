function v = observation_eq(Z,X,Y,a,gamma,D1,x,eta,delta,k1,k2,mu1,mu2)

b = stimuli_eq(D1,x,gamma,eta,delta);

v=sum(cell2mat(arrayfun(@(i)sum(squeeze(Z(i,:,:)).*(squeeze(log(twoPL(X(i,:),a,b))) - ones(size(Z,2),size(Z,3))*log(2*pi*expscale_besseli(k1)) + k1*cos(squeeze(Y(i,:,:))-ones(size(Z,2),size(Z,3))*mu1)) + ...
    (1-squeeze(Z(i,:,:))).*(squeeze(log(1-twoPL(X(i,:),a,b))) - ones(size(Z,2),size(Z,3))*log(2*pi*expscale_besseli(k2)) + k2*cos(squeeze(Y(i,:,:))-ones(size(Z,2),size(Z,3))*mu2)),1),1:size(X,1),'UniformOutput',false)));

end
