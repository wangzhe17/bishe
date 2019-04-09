clc;
clear all;
I=imread('56.png');
BW1=rgb2gray(I);           %BW1灰度图
thresh = graythresh(BW1);     %自动确定二值化阈值；
BW2 = im2bw(BW1,thresh);       %对图像自动二值化即可。 BW2二值图

BW2=imcomplement(BW2);    %颜色翻转
BW2 = bwareaopen(BW2,380,8);  %去文字
BW2=imcomplement(BW2);
figure,imshow(BW2);
BW3=BW2;
BW3=~BW3;

se=strel('disk',5);  
BW3=imclose(BW3,se);
BW3=~BW3;
A=imerode(BW3,se);     %一次腐蚀操作
[A1,AM]=zhang(A);
A=imerode(AM,se);     %一次腐蚀操作
[A1,AM]=zhang(A);
A=imerode(AM,se);     %一次腐蚀操作
[A1,AM]=zhang(A);
% A=imerode(AM,se);     %一次腐蚀操作
% [A1,AM]=zhang(A);
% A=imerode(AM,se);     %一次腐蚀操作
% [A1,AM]=zhang(A);
% A=imerode(AM,se);     %一次腐蚀操作
% [A1,AM]=zhang(A);
% 
% figure,imshow(AM);


A=AM;


se2=strel('disk',1);  
A=imerode(A,se2);     %一次腐蚀操作
edge1 = bwperim(A,8);       %edge1轮廓图\
% figure,imshow(AM);
% figure,imshow(A);
% figure,imshow(edge1);

[L,num] = bwlabel(edge1,8);
[w,h]=size(L);  %%%%L的宽高

images=regionprops(L,'Image');
status=regionprops(L,'BoundingBox');
centroid = regionprops(L,'Centroid');
%%
%%对每一个内部轮廓提取%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:num
    COPY=zeros(w,h);
    a1=ceil(status(i).BoundingBox(1));
    a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3));
    b1=ceil(status(i).BoundingBox(2));
    b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4));
    %%%%%%%%%%%%%%%单个连通域的提取%%%%%%%%%%%%%%%%%%%%%
    [r,c]=find(L==i);
    for x=1:size(r)
       COPY(r(x,1),c(x,1))=1; 
       B=COPY(b1:b2,a1:a2);
    end
    
    [r1,c1]=size(B);
    H=zeros(b2-b1+80,a2-a1+80);
    for j=1:r1
        for k=1:c1
            H(j+40,k+40)=B(j,k);
        end
    end
    xingz1{i}=H;
end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%对每一个单个连通域进行填充%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure();
for i=2:num
    X=xingz1{i};
    
    X=imfill(X,'holes');
    xingz2{i}=X;
%     subplot(7,7,i);
%     imshow(X);
end

%%%%%%%%检测凸包并填充%%%%%%%%%
% figure();
for i=2:num
    pic1=xingz1{i};
    [y2,x2]=find(pic1==1);
    dt = delaunayTriangulation(x2,y2);
    k = convexHull(dt);
    [r3,c3]=size(pic1);
    pic11=zeros(r3,c3);
    lkt = roipoly(pic11,x2(k),y2(k));%画出凸包的区域，取出来得到掩膜
    tuxingz{i}=lkt;
%    subplot(7,7,i);
%    imshow(lkt);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算凸包与原图形面积差%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:num
    picy=xingz2{i};
    pict=tuxingz{i};
    areay=regionprops(picy,'Area');
    areat=regionprops(pict,'Area');
    areay=areay.Area;
    areat=areat.Area;
    div=areat-areay;
    divs{i}=div;
end


% figure();
for i=2:num
    if divs{i}>100   %若面积差大于100
        continue;
    end
    houxuan=xingz2{i};
    houxuans{i}=houxuan;
%    subplot(7,7,i);
%    imshow(houxuan);
end
%%

%%
squares=zeros(1,num);
rectangles=zeros(1,num);
circles=zeros(1,num);
eclipses=zeros(1,num);
diamonds=zeros(1,num);
[cc,len]=size(houxuans);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%形状检测，采用图形面积、最小外接矩形面积+长宽比%%%%%%%%%%%%%%%%
for i=2:len
    houxuan=houxuans{i};
    if isempty(houxuan)
        continue;
    end
    [l,m] = bwlabel(houxuan);
    boundingBox=regionprops(l,'BoundingBox');
    area=regionprops(l,'Area');
    area=area.Area;
    a1=ceil(boundingBox.BoundingBox(1));
    a2=ceil(boundingBox.BoundingBox(1)+boundingBox.BoundingBox(3));
    b1=ceil(boundingBox.BoundingBox(2));
    b2=ceil(boundingBox.BoundingBox(2)+boundingBox.BoundingBox(4));
    areaR=(b2-b1)*(a2-a1);
    ratio=area/areaR;
    ckd=abs((b2-b1)-(a2-a1));%长宽比
    if(ratio>=0.95)
%         if(ckd<3)
%             squares(1,i)=1;
%         else
       rectangles(1,i)=1;   %只检测矩形
%         end
    end
    if(ratio>=0.73&&ratio<=0.83)
        if(ckd<6)
            circles(1,i)=1;
        else
            eclipses(1,i)=1;
        end
    end
    if(ratio>=0.45&&ratio<=0.55)
        diamonds(1,i)=1;
    end
    ratios{i}=ratio;
end
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%生成图元层%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w,h]=size(COPY); 
COPY2=ones(w,h);
for i=2:len
    houxuan=houxuans{i};
    if ~isempty(houxuan)
        a1=ceil(status(i).BoundingBox(1))+1;
        a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3));
        b1=ceil(status(i).BoundingBox(2))+1;
        b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4));
        image=images(i,1).Image;
        image=imcomplement(image);
        for j=b1:b2
            for z=a1:a2
                COPY2(j,z)=image(j-b1+1,z-a1+1); 
            end
        end
    end
end
% figure,imshow(COPY2);
% title('图元层');
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%生成连接线%%%%%%%%%%%%%%%%%%%%%%%%%%
COPY3=AM;
for i=2:len
    houxuan=houxuans{i};
    if ~isempty(houxuan)
        a1=ceil(status(i).BoundingBox(1))-2;
        a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3))+1;
        b1=ceil(status(i).BoundingBox(2))-2;
        b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4))+1;
        COPY3(b1:b2,a1:a2)=1;
    end
end
% figure,imshow(COPY3);
figure();
subplot(2,2,1),imshow(AM);
subplot(2,2,2),imshow(COPY2);
subplot(2,2,3),imshow(COPY3);
%%
COPY3=~COPY3;
[w2,h2]=size(COPY3); 
for i=2:w2-1
    for j=2:h2-1
        B=COPY3(i-1:i+1,j-1:j+1);
        if B(2,2)==1 && sum(B(:))==5
            if B(1,1)==1&&B(2,1)==1&&B(2,2)==1&&B(3,2)==1&&B(3,3)==1
                B(2,1)=0;
                B(3,2)=0;
                COPY3(i-1:i+1,j-1:j+1)=B;
            end
        end
    end
end
                    
figure,imshow(COPY3),title("test");
COPY3=~COPY3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
COPY3=imcomplement(COPY3);
[ll,mm] = bwlabel(COPY3,8);
figure,imshow(COPY3);

status2=regionprops(ll,'BoundingBox');
centroid2 = regionprops(ll,'Centroid');
i1=1;
title('连通域矩形切割');
for ts=1:mm
    ii=1;
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    for i1=b1:b2
        for j1=a1:a2
            B=COPY3(i1-1:i1+1,j1-1:j1+1);
            if B(2,2)==1                               %8邻域像素点唯一
                if sum(B(:))==2
                    K=[j1,i1];
                    Ks{ii}=K;
                    ii=ii+1;
                    endpoints{ts}=Ks;
                end
            end
        end
    end
    Ks=[];
end

for ts=1:mm
    ii=1;
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    for i1=b1:b2
        for j1=a1:a2
            B=COPY3(i1-1:i1+1,j1-1:j1+1);
            if B(2,2)==1                               %8邻域像素点为3
                if (B(2,1)==1&&B(1,2)==1&&B(2,3)==1)||(B(2,1)==1&&B(3,2)==1&&B(2,3)==1)||(B(1,2)==1&&B(2,1)==1&&B(3,2)==1)||(B(1,2)==1&&B(2,3)==1&&B(3,2)==1)||(B(1,2)==1&&B(2,3)==1&&B(3,1)==1)||(B(1,1)==1&&B(2,3)==1&&B(3,2)==1)
                    KK=[j1,i1];
                    KKs{ii}=KK;
                    ii=ii+1;
                    jointpoints{ts}=KKs;
                end
            end
        end
    end
    KKs=[];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,imshow(I);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%在原图中进行标注%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:num
    if(rectangles(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','c','LineWidth',3);
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '矩形','Color', 'c') 
    end
    if(diamonds(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','b');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '菱形','Color', 'b') 
    end
    if(circles(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','g','LineWidth',3);
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '圆形','Color', 'g') 
    end
    if(eclipses(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','c');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '椭圆','Color', 'c') 
    end
end
%%%%%%%%%%%%%%%画出端点位置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(endpoints)
    zuobiaos=endpoints{i};
    for j=1:length(zuobiaos)
        jj=zuobiaos{j};
        x0=jj(1,1);
        y0=jj(1,2);
        plot(x0,y0,'ro','LineWidth',5,'MarkerFaceColor','r')
    end
end
for i=1:length(jointpoints)
    zuobiaos=jointpoints{i};
    for j=1:length(zuobiaos)
        jj=zuobiaos{j};
        x0=jj(1,1);
        y0=jj(1,2);
        plot(x0,y0,'bo','LineWidth',5,'MarkerFaceColor','b')
    end
end
hold off;