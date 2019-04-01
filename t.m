clc;
clear all;
I=imread('xian.jpg');
BW1=I;
thresh = graythresh(BW1);     %�Զ�ȷ����ֵ����ֵ��
BW2 = im2bw(BW1,thresh);       %��ͼ���Զ���ֵ�����ɡ� BW2��ֵͼ

BW2=imcomplement(BW2);    %��ɫ��ת
[l,m] = bwlabel(BW2,8);
figure,imshow(BW2);

status=regionprops(l,'BoundingBox');
centroid = regionprops(l,'Centroid');
for i=1:m
    rectangle('position',status(i).BoundingBox,'edgecolor','r');
    text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), num2str(i),'Color', 'r') 
end
title('��ͨ������и�');

a1=ceil(status(2).BoundingBox(1));
a2=ceil(status(2).BoundingBox(1)+status(2).BoundingBox(3));
b1=ceil(status(2).BoundingBox(2));
b2=ceil(status(2).BoundingBox(2)+status(2).BoundingBox(4));

ii=1;
for i=b1:b2
    for j=a1:a2
        B=BW2(i-1:i+1,j-1:j+1);
        if B(2,2)==1;
            if sum(B(:))==2
                K=[j,i];
                Ks{ii}=K;
                ii=ii+1;
            end
        end
    end
end

