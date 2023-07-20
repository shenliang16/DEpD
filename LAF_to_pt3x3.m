function [Xp,Yp] = LAF_to_pt3x3(Xf, Yf)

if isempty(Xf)
   Xp = Xf; Yp=Yf; 
   return; 
end

if size(Xf,2)==6
    temp = sum(Xf);
    if temp(3)>0.05*temp(5) 
        Xp = Xf;
        Yp = Yf;  
        return;
    end
    affptx(:).x=Xf(:,5)';           affpty(:).x=Yf(:,5)';
    affptx(:).y=Xf(:,6)';           affpty(:).y=Yf(:,6)';
    affptx(:).a11=Xf(:,1)';         affpty(:).a11=Yf(:,1)';
    affptx(:).a12=Xf(:,3)';         affpty(:).a12=Yf(:,3)';
    affptx(:).a21=Xf(:,2)';         affpty(:).a21=Yf(:,2)';
    affptx(:).a22=Xf(:,4)';         affpty(:).a22=Yf(:,4)';

    X_3Pts=affpt_to_pt3x3(affptx);  Y_3Pts=affpt_to_pt3x3(affpty);
    % 将X和Y组合为一个矩阵用于输入：[Xpt1;Ypt1;Xpt2;Ypt2;Xpt3;Ypt3]
    Xp = X_3Pts([1:2,4:5,7:8],:)';   
    Yp = Y_3Pts([1:2,4:5,7:8],:)';  
    
else
%     temp = sum(Xf);
%     if temp(3)>0.05*temp(1)  && temp(4)>0.05*temp(1) && (sum(Xf(:,3)==31) + sum(Xf(:,3)-37.2<0.05))<5
%         Xp = Xf;
%         Yp = Yf;  
%         return;
%     end
    Xp1x=Xf(:,1);               Yp1x=Yf(:,1);
    Xp1y=Xf(:,2);               Yp1y=Yf(:,2);
    Sx=Xf(:,3);                 Sy=Yf(:,3);
    Rx=Xf(:,4);                 Ry=Yf(:,4);
    %% 角度是指逆时针方向，是否需要取反？    ORB
%     Rx=-Xf(:,4);                 Ry=-Yf(:,4);
    
%     temp = max(abs(Rx));
%     if temp>20
%         Rx = deg2rad(Rx);
%         Ry = deg2rad(Ry);  
%     end
   
%     if (sum(Xf(:,3)==31) + sum(abs(Xf(:,3)-37.2)<0.05))>5
%         Sx = Sx/2;
%         Sy = Sy/2;
%     end
    
    Xp2x = Xp1x + Sx.*sin(Rx); % sin：我们预先改变了X与Y的顺序，用X表示行，Y表示列 （VLfeat正好相反）
    Xp2y = Xp1y + Sx.*cos(Rx);
    Yp2x = Yp1x + Sy.*sin(Ry); % sin：我们预先改变了X与Y的顺序，用X表示行，Y表示列 （VLfeat正好相反）
    Yp2y = Yp1y + Sy.*cos(Ry); 
    
    Xp3x = Xp1x - Sx.*sin(Rx); % sin：我们预先改变了X与Y的顺序，用X表示行，Y表示列 （VLfeat正好相反）
    Xp3y = Xp1y - Sx.*cos(Rx);
    Yp3x = Yp1x - Sy.*sin(Ry); % sin：我们预先改变了X与Y的顺序，用X表示行，Y表示列 （VLfeat正好相反）
    Yp3y = Yp1y - Sy.*cos(Ry); 

%     Xp = [Xp1x, Xp1y, Xp2x, Xp2y, Xp3x, Xp3y];   
%     Yp = [Yp1x, Yp1y, Yp2x, Yp2y, Yp3x, Yp3y];  
    Xp = [Xp1x, Xp1y, Xp2x, Xp2y];   
    Yp = [Yp1x, Yp1y, Yp2x, Yp2y];
end


