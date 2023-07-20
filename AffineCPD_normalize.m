

% function  [X, Y, normal] =AffineCPD_normalize(x,y)
% global loc2
% if nargin<2, error('cpd_normalize error! Not enough input parameters.'); end;
% 
% if size(x,2)<=3
%     
%     n=size(x,1);
%     m=size(y,1);
% 
%     %% 
%     normal.xd=mean(x); normal.yd=mean(y); %
%     x=x-repmat(normal.xd,n,1);
%     y=y-repmat(normal.yd,m,1);
% 
%     %% 
%     normal.xscale=sqrt(sum(sum(x.^2,2))/n); 
%     normal.yscale=sqrt(sum(sum(y.^2,2))/m);
% 
%     X=x/normal.xscale;
%     Y=y/normal.yscale;
% 
% elseif size(x,2)==6
% 
%     n=size(x,1); m=size(y,1);
%     x1=x(:,1:2); x2=x(:,3:4); x3=x(:,5:6);
%     y1=y(:,1:2); y2=y(:,3:4); y3=y(:,5:6);
%     %% 
%     normal.x1d=mean(x1); normal.y1d=mean(y1); %
%     normal.x2d=mean(x2); normal.y2d=mean(y2); %
%     normal.x3d=mean(x3); normal.y3d=mean(y3); %
% 
%     x3=x3-repmat(normal.x3d,n,1);
%     y3=y3-repmat(normal.y3d,m,1);
%     % loc2(:,5:6)=loc2(:,5:6)-repmat(normal.y3d,m,1);
% 
%     normal.x1scale=sqrt(sum(sum(x1.^2,2))/n);
%     normal.x2scale=sqrt(sum(sum(x2.^2,2))/n);
%     normal.x3scale=sqrt(sum(sum(x3.^2,2))/n);
%     normal.y1scale=sqrt(sum(sum(y1.^2,2))/m);
%     normal.y2scale=sqrt(sum(sum(y2.^2,2))/m);
%     normal.y3scale=sqrt(sum(sum(y3.^2,2))/m);
% 
%     x1=x1/normal.x1scale;
%     x2=x2/normal.x2scale;
%     x3=x3/normal.x3scale;
%     y1=y1/normal.y1scale;
%     y2=y2/normal.y2scale;
%     y3=y3/normal.y3scale;
%     % loc2(:,1:4)=loc2(:,1:4)/normal.y1scale;
%     % loc2(:,5:6)=loc2(:,5:6)/normal.y2scale;
%     
% normal.FromA2T=[];
% %     normal.FromA2T=(normal.x1scale*normal.y3scale)/(normal.y1scale*normal.x3scale);
%     X=[x1,x2,x3]; Y=[y1,y2,y3];
% end




function  [X, Y, normal] =AffineCPD_normalize(x, y, flag, normal)
global loc2
if nargin<2, error('cpd_normalize error! Not enough input parameters.'); end;

if exist('normal', 'var')
    D = size(x, 2);  D_normal = size(normal.xd, 2);
    x = x - repmat(normal.xd, 1, D/D_normal);
    y = y - repmat(normal.yd, 1, D/D_normal);

    X=x./normal.xscale;
    Y=y./normal.yscale;

elseif exist('flag', 'var')
    %% reshape
    D = size(x,2);
    
    %% È¥µôNaN
    x_ = x';             x_(isnan(x_)) = [];
    y_ = y';             y_(isnan(y_)) = [];
    
    x_reshape = reshape(x_, 2, [])';
    y_reshape = reshape(y_, 2, [])';
    
    n=size(x_reshape,1);
    m=size(y_reshape,1);

    %% 
    normal.xd = mean(x_reshape);   
    normal.yd = mean(y_reshape); %
    
    x = x - repmat(normal.xd, 1, D/2);
    y = y - repmat(normal.yd, 1, D/2);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/m);

    X=x./normal.xscale;
    Y=y./normal.yscale;
    
elseif size(x,2)==2
    
    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean(x);   normal.yd=mean(y); %
    
    x=x-repmat(normal.xd,n,1);
    y=y-repmat(normal.yd,m,1);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;

elseif size(x,2)==4
    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean([x(:,1:2); x(:,3:4)]);   normal.yd=mean([y(:,1:2); y(:,3:4)]); %
    
    x=x-repmat(normal.xd,n,2);
    y=y-repmat(normal.yd,m,2);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/2/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/2/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;
    
elseif size(x,2)==6

    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean([x(:,1:2); x(:,3:4); x(:,5:6)]);   normal.yd=mean([y(:,1:2); y(:,3:4); y(:,5:6)]); %
    
    x=x-repmat(normal.xd,n,3);
    y=y-repmat(normal.yd,m,3);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/3/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/3/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;
    
    
    
%%   
% elseif size(x,2)==6
% 
%     n=size(x,1); m=size(y,1);
%     x1=x(:,1:4); x2=x(:,5:6);
%     y1=y(:,1:4); y2=y(:,5:6);
%     %% 
%     normal.x1d=mean(x1); normal.y1d=mean(y1); %
%     normal.x3d=mean(x2); normal.y3d=mean(y2); %
% 
%     x2=x2-repmat(normal.x3d,n,1);
%     y2=y2-repmat(normal.y3d,m,1);
%     % loc2(:,5:6)=loc2(:,5:6)-repmat(normal.y3d,m,1);
% 
%     normal.x1scale=sqrt(sum(sum(x1.^2,2))/(2*n));
%     normal.x2scale=sqrt(sum(sum(x2.^2,2))/n);
%     normal.y1scale=sqrt(sum(sum(y1.^2,2))/(2*m));
%     normal.y2scale=sqrt(sum(sum(y2.^2,2))/m);
% 
%     x1=x1/normal.x1scale;
%     x2=x2/normal.x2scale;
%     y1=y1/normal.y1scale;
%     y2=y2/normal.y2scale;
%     % loc2(:,1:4)=loc2(:,1:4)/normal.y1scale;
%     % loc2(:,5:6)=loc2(:,5:6)/normal.y2scale;
%     
% % normal.FromA2T=(normal.x2scale*normal.y1scale)/(normal.y2scale*normal.x1scale);
%     normal.FromA2T=(normal.x1scale*normal.y2scale)/(normal.y1scale*normal.x2scale);
%     X=[x1,x2]; Y=[y1,y2];
end

