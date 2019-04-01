clc;
clear all;
I=imread('xian.jpg');
BW1=I;
thresh = graythresh(BW1);     %自动确定二值化阈值；
BW2 = im2bw(BW1,thresh);       %对图像自动二值化即可。 BW2二值图

BW2=imcomplement(BW2);    %颜色翻转
[l,m] = bwlabel(BW2,8);
figure,imshow(BW2);

status=regionprops(l,'BoundingBox');
centroid = regionprops(l,'Centroid');
for i=1:m
    rectangle('position',status(i).BoundingBox,'edgecolor','r');
    text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), num2str(i),'Color', 'r') 
end
title('连通域矩形切割');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%找到所有端点%%%%%%%%%%%%%%%%%%%%%%%%%%
for ts=1:m
    ii=1;
    a1=ceil(status(ts).BoundingBox(1));
    a2=ceil(status(ts).BoundingBox(1)+status(ts).BoundingBox(3));
    b1=ceil(status(ts).BoundingBox(2));
    b2=ceil(status(ts).BoundingBox(2)+status(ts).BoundingBox(4));
    for i=b1:b2
        for j=a1:a2
            B=BW2(i-1:i+1,j-1:j+1);
            if B(2,2)==1                               %8邻域像素点唯一
                if sum(B(:))==2
                    K=[j,i];
                    Ks{ii}=K;
                    ii=ii+1;
                    endpoints{ts}=Ks;
                end
            end
        end
    end
    Ks=[];
end

for ts=1:m
    ii=1;
    a1=ceil(status(ts).BoundingBox(1));
    a2=ceil(status(ts).BoundingBox(1)+status(ts).BoundingBox(3));
    b1=ceil(status(ts).BoundingBox(2));
    b2=ceil(status(ts).BoundingBox(2)+status(ts).BoundingBox(4));
    for i=b1:b2
        for j=a1:a2
            B=BW2(i-1:i+1,j-1:j+1);
            if B(2,2)==1                               %8邻域像素点为3
                if (B(2,1)==1&&B(1,2)==1&&B(2,3)==1)||(B(2,1)==1&&B(3,2)==1&&B(2,3)==1)||(B(1,2)==1&&B(2,1)==1&&B(3,2)==1)||(B(1,2)==1&&B(2,3)==1&&B(3,2)==1)
                    KK=[j,i];
                    KKs{ii}=KK;
                    ii=ii+1;
                    jointpoints{ts}=KKs;
                end
            end
        end
    end
    KKs=[];
end

figure,imshow(BW2);
%%%%%%%%%%%%%%%画出端点位置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:length(endpoints)
%     zuobiaos=endpoints{i};
%     for j=1:length(zuobiaos)
%         jj=zuobiaos{j};
%         x0=jj(1,1);
%         y0=jj(1,2);
%         plot(x0,y0,'ro')
%     end
% end