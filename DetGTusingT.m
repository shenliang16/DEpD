function [GrdTrth, Dis0, crsp, Dis]=DetGTusingT(X,Y,Im1,Im2,thr,T)

if nargin<6
    T=GetTransformation(Im1,Im2);
end
    N=size(X,1);    X3=[X(:,1:2),ones(N,1)];     
    XT=(T*X3')';        TransformedX=XT(:,1:2)./XT(:,3);
    
    Y = Y(1:N,:);
    Dis0 = Y(:,1:2)-TransformedX;
    Dis=sum(Dis0.^2, 2);
    GrdTrth=(sqrt(Dis)<thr);

    
    Dis1 = Y(:,1)-TransformedX(:,1)';
    Dis2 = Y(:,2)-TransformedX(:,2)';
    Distance = Dis1.^2 + Dis2.^2;
    [disInd, ind] = min(Distance,[],2);
    indGT = find(disInd<thr^2);
    crsp = [ind(indGT),indGT];
% sum(GrdTrth)