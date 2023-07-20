%% 单张图片

function Matching_Plot_PLR(Im11,Im22,Corresp,Flag,fileNm,loc1_input,loc2_input)
global loc1 loc2 X Y X_ind Y2_ind

if exist('loc1_input','var') && exist('loc2_input','var')
    loc1_Crrt=loc1_input;
    loc2_Crrt=loc2_input;
else
    loc1_Crrt=loc1;
    loc2_Crrt=loc2;
end

% if size(Im11,3)==3
%     Im11 = rgb2gray(Im11);
% end
% if size(Im22,3)==3
%     Im22 = rgb2gray(Im22);
% end
if size(Im11,3)==3 && size(Im22,3)==1
    Im11 = rgb2gray(Im11);
end
if size(Im22,3)==3 && size(Im11,3)==1
    Im22 = rgb2gray(Im22);
end

% 统一大小
[m1,n1,~]=size(Im11);  [m2,n2,~]=size(Im22);  sizeMin=max([m1,n1],[m2,n2]);  
size3=size(Im11,3);
Im1=zeros([sizeMin,size3]);
Im2=zeros([sizeMin,size3]);

Im1(1:m1,1:n1,:)=Im11; 
Im2(1:m2,1:n2,:)=Im22; 


if max(Im1(1))<=1
    Im1=uint8(255*Im1);
    Im2=uint8(255*Im2);
end
interval = 20; 
WhiteInterval = uint8(255*ones(size(Im1,1), interval, size3));
figure('position',[0 0 1000 400]); %set(gcf,'position',[0 0 500 500])
imagesc(cat(2, Im1, WhiteInterval, Im2)); 
colormap('gray'); 
hold on; 
% cols1 = 3*size(Im1,2)/5;
cols1 = size(Im1,2);
% times=[1,1,5,5,5,5];
times=1;

plusNum=[cols1+interval,0,0,0,0,0];

temp = sum(loc1_Crrt(Corresp(:,1),:),1);

if isempty(Corresp);
    axis equal; axis off;
    return;
end


if temp(3)>0.05*temp(1)  && temp(4)>0.05*temp(1)
    TwoPts_flag = 1;
    Line1_2tps = loc1_Crrt(Corresp(:,1),1:4);
    Line2_2tps = loc2_Crrt(Corresp(:,2),1:4);
else
    TwoPts_flag = 0;
    [Line1_2tps, Line2_2tps] = LAF_to_pt3x3(loc1_Crrt(Corresp(:,1),1:4), loc2_Crrt(Corresp(:,2),1:4));
end

for i=1:size(Corresp,1) 
    data = loc1_Crrt(Corresp(i,1),:);
    D = length(data(~isnan(data)));   ind = D/2;
%     if all(~isnan(loc1_Crrt(Corresp(i,1),:))); times=[1,1,1.5,1.5,1.5,1.5]; end;
    if Flag(i)==-1
        color='w';
        if ~TwoPts_flag
            vl_plotframe_oursNoCircle(loc1_Crrt(Corresp(i,1),[2,1,6,5,4,3]).*times,color(ind,:)) ; 
            vl_plotframe_oursNoCircle((loc2_Crrt(Corresp(i,2),[2,1,6,5,4,3])+plusNum).*times,color(ind,:)) ; 
        else
            line([loc1_Crrt(Corresp(i,1),4)' ;  loc1_Crrt(Corresp(i,1),2)'], [loc1_Crrt(Corresp(i,1),3)'; loc1_Crrt(Corresp(i,1),1)'],'Color',color(ind,:),'LineWidth',2.5);
            line([loc2_Crrt(Corresp(i,2),4)';  loc2_Crrt(Corresp(i,2),2)']+cols1+interval, [loc2_Crrt(Corresp(i,2),3)'; loc2_Crrt(Corresp(i,2),1)'],'Color',color(ind,:),'LineWidth',2.5);
        end
        line([loc1_Crrt(Corresp(i,1),2)'; loc2_Crrt(Corresp(i,2),2)'+cols1+interval], [loc1_Crrt(Corresp(i,1),1)' ; loc2_Crrt(Corresp(i,2),1)'],'linestyle', ':','linewidth', 0.5, 'color', color,'Marker','.','Markersize',9) ;%'g'
    elseif Flag(i)==0
        color=[1,0,0; 0,1,1; 1,0,0];
%         Color_rnd_sel = rand(1,3);
%         color=[1,0,0; Color_rnd_sel; 0,0,1];
        if ~TwoPts_flag
            vl_plotframe_oursNoCircle(loc1_Crrt(Corresp(i,1),[2,1,6,5,4,3]).*times,color(ind,:)) ; 
            vl_plotframe_oursNoCircle((loc2_Crrt(Corresp(i,2),[2,1,6,5,4,3])+plusNum).*times,color(ind,:)) ; 
        else
            line([loc1_Crrt(Corresp(i,1),4)' ;  loc1_Crrt(Corresp(i,1),2)'], [loc1_Crrt(Corresp(i,1),3)'; loc1_Crrt(Corresp(i,1),1)'],'Color',color(ind,:),'LineWidth',2.5);
            line([loc2_Crrt(Corresp(i,2),4)';  loc2_Crrt(Corresp(i,2),2)']+cols1+interval, [loc2_Crrt(Corresp(i,2),3)'; loc2_Crrt(Corresp(i,2),1)'],'Color',color(ind,:),'LineWidth',2.5);
        end
%         color=[1,0,0; 1,0.0,0.0; 0,0,1];
        color=[1,0,0; 0,1,1; 0,0,1];
        line([loc1_Crrt(Corresp(i,1),2)'; loc2_Crrt(Corresp(i,2),2)'+cols1+interval], [loc1_Crrt(Corresp(i,1),1)' ;  loc2_Crrt(Corresp(i,2),1)'],'linestyle', '-','linewidth', 0.3, 'color',color(ind,:),'Marker','.','Markersize',9) ;%  [0.8,0.1,0],'Marker','.'
    elseif Flag(i)==1
%         Color_rnd_sel = rand(1,3);
%         color=[1,0,0; Color_rnd_sel; 0,0,1];
        color=[0,1,0; 1,0,0; 0,0,1];
        if ~TwoPts_flag
            vl_plotframe_oursNoCircle(loc1_Crrt(Corresp(i,1),[2,1,6,5,4,3]).*times,color(ind,:)) ; 
            vl_plotframe_oursNoCircle((loc2_Crrt(Corresp(i,2),[2,1,6,5,4,3])+plusNum).*times,color(ind,:)) ; 
        else
            line([loc1_Crrt(Corresp(i,1),4)' ;  loc1_Crrt(Corresp(i,1),2)'], [loc1_Crrt(Corresp(i,1),3)'; loc1_Crrt(Corresp(i,1),1)'],'Color',color(ind,:),'LineWidth',2.5);
            line([loc2_Crrt(Corresp(i,2),4)';  loc2_Crrt(Corresp(i,2),2)']+cols1+interval, [loc2_Crrt(Corresp(i,2),3)'; loc2_Crrt(Corresp(i,2),1)'],'Color',color(ind,:),'LineWidth',2.5);
        end
        % 添加连接线
        color=[1,0,0; 1,0,0; 0,0,1];
        line([loc1_Crrt(Corresp(i,1),2)'; loc2_Crrt(Corresp(i,2),2)'+cols1+interval], [loc1_Crrt(Corresp(i,1),1)' ;  loc2_Crrt(Corresp(i,2),1)'],'linestyle', '-','linewidth', 0.3, 'color',color(ind,:) ,'Marker','.','Markersize',9) ;%[0,0.5,0.8]%'Marker','.'
    end

% 显示编号
if Flag(i)>=0 %&& Flag(i)==0
    FontSize = 12;
    Mid_Loc = mean(reshape(Line1_2tps(i,[2,1,4,3]),2,2)');      Mid_Loc = double(Mid_Loc);
    text(Mid_Loc(1),Mid_Loc(2),num2str(i),'Color',color(ind,:),'FontSize',FontSize)
    Mid_Loc = mean(reshape(Line2_2tps(i,[2,1,4,3])+plusNum([1,2,1,2]),2,2)');       Mid_Loc = double(Mid_Loc);
    text(Mid_Loc(1),Mid_Loc(2),num2str(i),'Color',color(ind,:),'FontSize',FontSize)
end
end
 axis equal; axis off;
% 保存
if exist(fileNm,'var') && ~isempty(fileNm)
    axis equal; axis off;  cd ./Result;  saveas(gcf,fileNm); cd ..;
end
set(gca, 'LooseInset', [0,0,0,0]);
end
