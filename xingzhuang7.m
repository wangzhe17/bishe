clc;
clear all;
I=imread('56.png');
I=imresize(I,2);
BW1=rgb2gray(I);           %BW1灰度图
thresh = graythresh(BW1);     %自动确定二值化阈值；
BW2 = im2bw(BW1,thresh);       %对图像自动二值化即可。 BW2二值图

BW2=imcomplement(BW2);    %颜色翻转
BW2 = bwareaopen(BW2,1580,8);  %去文字
BW2=imcomplement(BW2);

BW3=BW2;
BW3=~BW3;

se=strel('disk',10);  
BW3=imclose(BW3,se);
BW3=~BW3;
% figure,imshow(BW3);
A=imerode(BW3,se);     %一次腐蚀操作
[A1,AM]=zhang(A);




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
for iz=2:len
    houxuan=houxuans{iz};
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
    if(ratio>=0.93)
%         if(ckd<3)
%             squares(1,i)=1;
%         else
       rectangles(1,iz)=1;   %只检测矩形
%         end
    end
    if(ratio>=0.66&&ratio<=0.87)
        if(ckd<10)
            circles(1,iz)=1;
        else
            eclipses(1,iz)=1;
        end
    end
    if(ratio>=0.45&&ratio<=0.55)
        diamonds(1,iz)=1;
    end
    ratios{iz}=ratio;
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
        if rectangles(1,i)==1
            a1=ceil(status(i).BoundingBox(1))-2;
            a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3))+1;
            b1=ceil(status(i).BoundingBox(2))-2;
            b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4))+1;
            COPY3(b1:b2,a1:a2)=1;
        end
    end
end

for i=2:len
    houxuan=houxuans{i};
    if ~isempty(houxuan)
        if circles(1,i)==1
            a1=ceil(status(i).BoundingBox(1))-4;
            a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3))+4;
            b1=ceil(status(i).BoundingBox(2))-4;
            b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4))+4;
            tmp=~COPY3(b1:b2,a1:a2);
            c1=(b2-b1)/2+1;
            c2=(a2-a1)/2+1;
            r=min([(b2-b1),(a2-a1)])/2;
            theta=0:2*pi/3600:2*pi;
            Circle1=c1+r*cos(theta);
            Circle2=c2+r*sin(theta);
            l=roipoly(tmp,Circle1,Circle2);
            for i1= 1:(b2-b1)
                for j1=1:(a2-a1)
                    if l(i1,j1)==1
                        tmp(i1,j1)=0;
                    end
                end
            end
            COPY3(b1:b2,a1:a2)=~tmp;
        end
    end
end
COPY3=imcomplement(COPY3);    %颜色翻转
COPY3 = bwareaopen(COPY3,50,8);  
COPY3=imcomplement(COPY3);
% figure,imshow(COPY3),title('连接线');
figure();
subplot(2,2,1),imshow(AM);
subplot(2,2,2),imshow(COPY2);
subplot(2,2,3),imshow(COPY3);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%连接线箭头
COPY5=BW2;
for i=2:len
    houxuan=houxuans{i};
    if ~isempty(houxuan)
        if rectangles(1,i)==1
            a1=ceil(status(i).BoundingBox(1))-6;
            a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3))+7;
            b1=ceil(status(i).BoundingBox(2))-6;
            b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4))+7;
            COPY5(b1:b2,a1:a2)=1;
        end
    end
end

for i=2:len
    houxuan=houxuans{i};
    if ~isempty(houxuan)
        if circles(1,i)==1
            a1=ceil(status(i).BoundingBox(1))-9;
            a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3))+9;
            b1=ceil(status(i).BoundingBox(2))-9;
            b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4))+9;
            tmp=~COPY5(b1:b2,a1:a2);
            c1=(b2-b1)/2+1;
            c2=(a2-a1)/2+1;
            r=min([(b2-b1),(a2-a1)])/2;
            theta=0:2*pi/3600:2*pi;
            Circle1=c1+r*cos(theta);
            Circle2=c2+r*sin(theta);
            l=roipoly(tmp,Circle1,Circle2);
%             figure,imshow(l);
            for i1= 1:(b2-b1)
                for j1=1:(a2-a1)
                    if l(i1,j1)==1
                        tmp(i1,j1)=0;
                    end
                end
            end
            COPY5(b1:b2,a1:a2)=~tmp;
%             figure,imshow(tmp);
        end
    end
end
% lkt = roipoly(pic11,x2(k),y2(k))
COPY5=imcomplement(COPY5);    %颜色翻转
COPY5 = bwareaopen(COPY5,180,8);  
COPY5=imcomplement(COPY5);
figure,imshow(COPY5),title('连接线');

se1=strel('disk',10);  
COPY6=~COPY5;
COPY6=imclose(COPY6,se1);
COPY6=~COPY6;
AA=imerode(COPY6,se1);     
[AA1,AAM]=zhang(AA);
AA=imerode(AAM,se1);     
[AA1,AAM]=zhang(AA);
figure,imshow(AAM)




%%
COPY6=AAM;
%%
%%进一步细化
COPY6=~COPY6;
[w2,h2]=size(COPY6); 
for i=2:w2-1
    for j=2:h2-1
        B=COPY6(i-1:i+1,j-1:j+1);
        if B(2,2)==1 && sum(B(:))==5
            if B(1,1)==1&&B(2,1)==1&&B(2,2)==1&&B(3,2)==1&&B(3,3)==1
                B(2,1)=0;
                B(3,2)=0;
                COPY6(i-1:i+1,j-1:j+1)=B;
            end
        end
        if B(2,2)==1 && sum(B(:))==5
            if B(1,3)==1&&B(2,2)==1&&B(2,3)==1&&B(3,1)==1&&B(3,2)==1
                B(2,3)=0;
                B(3,2)=0;
                COPY6(i-1:i+1,j-1:j+1)=B;
            end
        end
        if B(2,2)==1 && sum(B(:))==4
            if B(1,1)==1&&B(2,1)==1&&B(2,2)==1&&B(3,3)==1
                B(2,1)=0;
                COPY6(i-1:i+1,j-1:j+1)=B;
            end
        end
        if B(2,2)==1 && sum(B(:))==4
            if B(1,3)==1&&B(2,2)==1&&B(2,3)==1&&B(3,1)==1
                B(2,3)=0;
                COPY6(i-1:i+1,j-1:j+1)=B;
            end
        end
    end
end                   
% figure,imshow(COPY3),title("test");
COPY6=~COPY6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
COPY6=imcomplement(COPY6);
[ll,mm] = bwlabel(COPY6,8);
figure,imshow(COPY6);

status2=regionprops(ll,'BoundingBox');
centroid2 = regionprops(ll,'Centroid');
i1=1;

for ts=1:mm
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    A=zeros(w,h);
    [r,c]=find(ll==ts);
    for x=1:size(r)
        A(r(x,1),c(x,1))=1;
    end
    endpointlines{ts}=A;
end

for ts=1:mm
    ii=1;
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    A=endpointlines{ts};
    for i1=b1:b2
        for j1=a1:a2
            B=A(i1-1:i1+1,j1-1:j1+1);
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
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    A=zeros(w,h);
    [r,c]=find(ll==ts);
    for x=1:size(r)
        A(r(x,1),c(x,1))=1;
    end
    jointpointlines{ts}=A;
end
for ts=1:mm
    ii=1;
    a1=ceil(status2(ts).BoundingBox(1));
    a2=ceil(status2(ts).BoundingBox(1)+status2(ts).BoundingBox(3));
    b1=ceil(status2(ts).BoundingBox(2));
    b2=ceil(status2(ts).BoundingBox(2)+status2(ts).BoundingBox(4));
    A=jointpointlines{ts};
    for i1=b1:b2
        for j1=a1:a2
            B=A(i1-1:i1+1,j1-1:j1+1);
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
J=~COPY5;


[L,num1] = bwlabel(J,8);
[w,h]=size(J);  %%%%L的宽高
lines=regionprops(J,'BoundingBox');


for i=1:num1
    T=zeros(w,h);
    a1=ceil(lines(i).BoundingBox(1));
    a2=ceil(lines(i).BoundingBox(1)+lines(i).BoundingBox(3));
    b1=ceil(lines(i).BoundingBox(2));
    b2=ceil(lines(i).BoundingBox(2)+lines(i).BoundingBox(4));
    A=zeros(b2-b1,a2-a1);
    [r,c]=find(L==i);
    for x=1:size(r)
        T(r(x,1),c(x,1))=1;
    end
%     figure,imshow(T);
    arrowlines{i}=T;
end


%%


%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,imshow(I);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%在原图中进行标注%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ia=1:num
    if(rectangles(1,ia)==1)
        rectangle('position',status(ia).BoundingBox,'edgecolor','c','LineWidth',5);
        text(centroid(ia,1).Centroid(1,1),centroid(ia,1).Centroid(1,2), '矩形','Color', 'c') 
    end
    if(diamonds(1,ia)==1)
        rectangle('position',status(ia).BoundingBox,'edgecolor','b');
        text(centroid(ia,1).Centroid(1,1),centroid(ia,1).Centroid(1,2), '菱形','Color', 'b') 
    end
    if(circles(1,ia)==1)
        rectangle('position',status(ia).BoundingBox,'Curvature',[1,1],'edgecolor','g','LineWidth',5);
        text(centroid(ia,1).Centroid(1,1),centroid(ia,1).Centroid(1,2), '圆形','Color', 'g') 
    end
    if(eclipses(1,ia)==1)
        rectangle('position',status(ia).BoundingBox,'edgecolor','c');
        text(centroid(ia,1).Centroid(1,1),centroid(ia,1).Centroid(1,2), '椭圆','Color', 'c') 
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
if(~isempty(jointpoints))
    for i=1:length(jointpoints)
        zuobiaos=jointpoints{i};
        for j=1:length(zuobiaos)
            jj=zuobiaos{j};
            x0=jj(1,1);
            y0=jj(1,2);
            plot(x0,y0,'bo','LineWidth',5,'MarkerFaceColor','b')
        end
    end
end
hold off;

%%

