I=imread('20.jpg');
BW1=rgb2gray(I);           %BW1�Ҷ�ͼ
thresh = graythresh(BW1);     %�Զ�ȷ����ֵ����ֵ��
BW2 = im2bw(BW1,thresh);       %��ͼ���Զ���ֵ�����ɡ� BW2��ֵͼ

BW2=imcomplement(BW2);    %��ɫ��ת
BW2 = bwareaopen(BW2,200,8);  %ȥ����
BW2=imcomplement(BW2);

se=strel('disk',1);  
A=imerode(BW2,se);     %һ�θ�ʴ����
A1=imdilate(A,se);
edge1 = bwperim(A,8);       %edge1����ͼ

[L,num] = bwlabel(edge1,8);

[w,h]=size(L);  %%%%L�Ŀ��

status=regionprops(L,'BoundingBox');
centroid = regionprops(L,'Centroid');
%%%%%%%%%%��ÿһ���ڲ�������ȡ%%%%%%%%%%%%%%%%%
for i=2:num
    COPY=zeros(w,h);
    a1=ceil(status(i).BoundingBox(1));
    a2=ceil(status(i).BoundingBox(1)+status(i).BoundingBox(3));
    b1=ceil(status(i).BoundingBox(2));
    b2=ceil(status(i).BoundingBox(2)+status(i).BoundingBox(4));
    %%%%%%%%%%%%%%%������ͨ�����ȡ%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%��ÿһ��������ͨ��������%%%%%%
% figure();
for i=2:num
    X=xingz1{i};
    
    X=imfill(X,'holes');
    xingz2{i}=X;
%     subplot(7,7,i);
%     imshow(X);
end

%%%%%%%%���͹�������%%%%%%%%%
% figure();
for i=2:num
    pic1=xingz1{i};
    [y2,x2]=find(pic1==1);
    dt = delaunayTriangulation(x2,y2);
    k = convexHull(dt);
    [r3,c3]=size(pic1);
    pic11=zeros(r3,c3);
    lkt = roipoly(pic11,x2(k),y2(k));%����͹��������ȡ�����õ���Ĥ
    tuxingz{i}=lkt;
%     subplot(7,7,i);
%     imshow(lkt);
end

%%%%%%%%����͹����ԭͼ�������%%%%%%%%%%%
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
    if divs{i}>0   %����������0
        continue;
    end
    houxuan=xingz2{i};
    houxuans{i}=houxuan;
%     subplot(7,7,i);
%     imshow(houxuan);
end



squares=zeros(1,num);
rectangles=zeros(1,num);
circles=zeros(1,num);
eclipses=zeros(1,num);
diamonds=zeros(1,num);
[cc,len]=size(houxuans);
%%%%%%%%%%%%%%��״��⣬����ͼ���������С��Ӿ������+�����%%%%%%%%%%%%%%%%
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
    ckd=abs((b2-b1)-(a2-a1));%�����
    if(ratio>=0.95)
        if(ckd<3)
            squares(1,i)=1;
        else
            rectangles(1,i)=1;
        end
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

%%%%%%%%%%%%%%��ԭͼ�н��б�ע%%%%%%%%%%%%%%%%
figure,imshow(BW2)
for i=2:num
    if(rectangles(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','r');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '����','Color', 'r') 
    end
    if(diamonds(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','b');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '����','Color', 'b') 
    end
    if(circles(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','g');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), 'Բ��','Color', 'g') 
    end
    if(eclipses(1,i)==1)
        rectangle('position',status(i).BoundingBox,'edgecolor','c');
        text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), '��Բ','Color', 'c') 
    end
end
