function [A,M]=zhang(filename)
tic();
% A=imread (filename) ; 
% A=im2bw (A,0.8);
A=filename;
I=~A ; 
[m n]=size (I) ; 
%Zhang Quick Method
for nc=1:60
for i=1:m
    for j=1:n 
        if I(i,j)==0
            M(i,j)=0; 
        else 
            P(1)=I(i,j) ; 
            if i-1>0 
                P(2) = I(i-1,j) ;
            else
                P(2) = 0 ; 
            end
            if (j+1)<=n & (i-1)>0 
                P(3) = I(i-1,j+1) ;
            else
                P(3) = 0 ; 
            end
            if (j+1)<=n 
                P(4) = I(i,j+1) ;
            else
                P(4) = 0 ; 
            end
            if (i+1)<=m & (j+1)<=n 
                P(5) = I(i+1,j+1) ;
            else
                P(5) = 0 ; 
            end
            if (i+1)<=m
                P(6)= I(i+1,j) ;
            else
                P(6) = 0 ; 
            end
            if (j-1)>0 & (i+1)<=m
                P(7) = I(i+1,j-1) ;
            else
                P(7) = 0 ; 
            end
            if (j-1)>0 
                P(8) = I(i,j-1) ;
            else
                P(8) = 0 ; 
            end
            if (i-1)>0 & (j-1)>0 
                P(9) = I(i-1,j-1) ;
            else
                P(9) = 0 ; 
            end
            P(10)=P(1) ; 
            Np1=sum(P(2:9)) ;
            Sp1=0; 
            Sp1=abs(P(9)-P(2)) ; 
            for t=2:8
                 Sp1=Sp1+abs(P(t)-P(t+1));
            end
            Sp1=Sp1/2 ; 
            if mod(nc,2)~=0 
             if (Np1>=2) & (Np1<=6) & (Sp1==1) & (P(2)*P(4)*P(6)==0) & (P(4)*P(6)*P(8)==0) 
                   M(i,j) = 0 ;
              else 
                   M(i,j) = 1 ;
               end
            else
             if Np1>=2 & Np1<=6 & Sp1==1 & P(2)*P(4)*P(8)==0 & P(2)*P(6)*P(8)==0 
                   M(i,j) = 0 ;
              else 
                   M(i,j) = 1 ;
              end
    
            end
            
        end
    end
end
I=M; 
end
M=~M;
timespent=toc
end