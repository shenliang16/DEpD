function [Np_pts,Np_line,P1_pts,P1_line,E_pts,E_line]=AffineCPD_Pts_Line(  Xi2,Ti2,  Xi4,Ti4, sigma2, omiga)


[N2, D2]=size(Xi2);       [N4, D4]=size(Xi4);     u = 1;    
    
%% 2D point
    if isempty(Xi2)
        P1_2 = [];  E2 = 0;    if length(omiga)==1;    omiga = repmat(omiga, 1,2); end
    else
        X1=Xi2(:,1:2);
        T1=Ti2(:,1:2);

        son2 = ((2*pi)^(D2/2)) * (sigma2^(D2/2)) * (omiga(1)/(1-omiga(1))) / u^(D2/2);

        sontemp = (X1-T1).*(X1-T1)/(2*sigma2); 
        sontemp1 = sum(sontemp,2); 
        son1 = exp(-sontemp1);  % 概率向量
        son = son1+son2;
        P1_2 = son1./son;
        E2 = sum(P1_2.*sontemp1);

        Np=sum(P1_2);  
        E2 = E2 + Np*log(sigma2);
        E2 = E2 - Np*log((1-omiga(1)))-(N2-Np)*log(omiga(1));  
    end
%% 4D line
    if isempty(Xi4)
        P1_4 = [];  E4 = [];
    else
        X1=Xi4(:,1:2);   X2=Xi4(:,3:4);
        T1=Ti4(:,1:2);   T2=Ti4(:,3:4);

        son2 = ((2*pi)^(D4/2)) * (sigma2^(D4/2)) * (omiga(2)/(1-omiga(2))) / u^(D4/2);

        sontemp = (X1-T1).*(X1-T1)/(2*sigma2)+...
                  (X2-T2).*(X2-T2)/(2*sigma2); 
        sontemp1 = sum(sontemp,2); 
        son1 = exp(-sontemp1);  % 概率向量
        son = son1+son2;
        P1_4 = son1./son;
        E4 = sum(P1_4.*sontemp1);

        Np=sum(P1_4);  
        E4 = E4 + 2*Np*log(sigma2);
        E4 = E4 - Np*log((1-omiga(2)))-(N4-Np)*log(omiga(2));   
    end

    %% Combine    
    Np_pts = sum([P1_2]);      
    Np_line = sum([P1_4]);
    P1_pts = diag(P1_2);
    P1_line = diag(P1_4);
    E_pts = E2;
    E_line = E4;
end