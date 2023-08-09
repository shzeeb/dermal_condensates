%triangualr lattice creator: generates a triangualr latttice of hair
%follicles with hexagonal packing. Inputs:
%Lx and Ly: domain sizes in x and y directions
%num_fols: number of follicles (nodes) on the lattice
%r: radius of a follicle, same for all follicles
%d_sep: distace which separated centre of one follicle form that of
%adjacent one.

%Created by: Shahzeb Raja Noureen
%Date created: 13/12/2020
%Last updated: 13/12/2020


function [xfc,yfc]=follicles(Lx,Ly,num_fol,r,d_sep,ghost)

xfc=zeros(num_fol,1);
yfc=zeros(num_fol,1);

%set d_x and d_y to 0 for mistaligned lattice (not centrered)
d_x=30;
d_y=60;

%follicle arrangement and their centre coordinates for different number of
%follicles
if num_fol==1
    xfc(1)=Lx/2;
    yfc(1)=Ly/2;
    
elseif num_fol==4
    
    xfc(1) = 158.5;
    xfc(2) = 158.5;
    xfc(3) = 475.5;
    xfc(4) = 475.5;

    yfc(1) = 158.5;
    yfc(2) = 475.5;
    yfc(3) = 475.5;
    yfc(4) = 158.5;


%coordiantes if num_fol=9
elseif num_fol==9
    
    xfc(1) = 170.75;
    xfc(2) = xfc(1) + d_sep;
    xfc(3) = xfc(2) + d_sep;
    xfc(4) = 170.75+84.56;
    xfc(5) = xfc(4) + d_sep;
    xfc(6) = xfc(5) + d_sep;
    xfc(7) = xfc(1);
    xfc(8) = xfc(2);
    xfc(9) = xfc(3);


    yfc(1) = 170.95;
    yfc(2) = 170.95;
    yfc(3) = 170.95;
    yfc(4) = 317.08;
    yfc(5) = 317.08;
    yfc(6) = 317.08;
    yfc(7) = 317.08+146.15;
    yfc(8) = 317.08+146.15;
    yfc(9) = 317.08+146.15;

%coordiantes of follicles centres if num_fol=14
elseif num_fol==14

    xfc(1)=r+d_x;
    xfc(2)=xfc(1)+d_sep;
    xfc(3)=xfc(1)+2*d_sep;
    xfc(4)=xfc(1)+3*d_sep;

    xfc(5)=r+d_sep/2+d_x;
    xfc(6)=xfc(5)+d_sep;
    xfc(7)=xfc(5)+2*d_sep;

    xfc(8)=r+d_x;
    xfc(9)=xfc(1)+d_sep;
    xfc(10)=xfc(1)+2*d_sep;
    xfc(11)=xfc(1)+3*d_sep;

    xfc(12)=r+d_sep/2+d_x;
    xfc(13)=xfc(5)+d_sep;
    xfc(14)=xfc(5)+2*d_sep;




    yfc(1:4)=r+d_y;


    yfc(5:7)=yfc(1)+sqrt(170*170-85*85);%d_sep;


    yfc(8:11)=yfc(1)+2*sqrt(170*170-85*85);%d_sep;


    yfc(12:14)=yfc(1)+3*sqrt(170*170-85*85);%d_sep;

    
    %ghose follicles: incorrect interpretaton: redo. just reflect the main
    %follicles in and x and y axes
    if ghost==1
        
        %adjacent x centres
        xfc([15 16 17 18])=[xfc(1) xfc(5) xfc(8) xfc(12)]+Lx;
        
        xfc([19 20 21 22])=[xfc(4) xfc(7) xfc(11) xfc(14)]-Lx;
        
        xfc([23 24 25 26])=[xfc(1) xfc(2) xfc(3) xfc(4)];
        
        xfc([27 28 29])=[xfc(12) xfc(13) xfc(14)];
        
        %diagonal x centres
        xfc([30 31 32 33])=[xfc(1)+Lx xfc(4)-Lx xfc(14)-Lx xfc(12)+Lx];
        
        
        %adjacent y centres
        yfc([15 16 17 18])=[yfc(1) yfc(5) yfc(8) yfc(12)];
        
        yfc([19 20 21 22])=[yfc(4) yfc(7) yfc(11) yfc(14)];
        
        yfc([23 24 25 26])=[yfc(1) yfc(2) yfc(3) yfc(4)]+Ly;
        
        yfc([27 28 29])=[yfc(12) yfc(13) yfc(14)]-Ly;
        
        %diagonal y centres
        yfc([30 31 32 33])=[yfc(1)+Ly yfc(4)+Ly yfc(14)-Ly yfc(12)-Ly];    
        
    end
    

%coordiantes of follicles centres if num_fol=16
elseif num_fol==16

    xfc(1)=r;
    xfc(2)=xfc(1)+d_sep;
    xfc(3)=xfc(1)+2*d_sep;
    xfc(4)=xfc(1)+3*d_sep;

    xfc(5)=r+d_sep/2;
    xfc(6)=xfc(5)+d_sep;
    xfc(7)=xfc(5)+2*d_sep;
    xfc(8)=xfc(5)+3*d_sep;

    xfc(9)=r;
    xfc(10)=xfc(1)+d_sep;
    xfc(11)=xfc(1)+2*d_sep;
    xfc(12)=xfc(1)+3*d_sep;

    xfc(13)=r+d_sep/2;
    xfc(14)=xfc(5)+d_sep;
    xfc(15)=xfc(5)+2*d_sep;
    xfc(16)=xfc(5)+3*d_sep;



    yfc(1:4)=r+60;

    yfc(5:8)=yfc(1)+sqrt(170*170-85*85);%d_sep;

    yfc(9:12)=yfc(1)+2*sqrt(170*170-85*85);%d_sep;

    yfc(13:16)=yfc(1)+3*sqrt(170*170-85*85);%d_sep;

    
    %ghost nodes
    if ghost==1
        xfc(17:32)=xfc(1:16)+Lx;

        xfc(33:48)=xfc(1:16);

        xfc(49:64)=xfc(1:16)-Lx;

        xfc(65:80)=xfc(1:16);


        yfc(17:32)=yfc(1:16);

        yfc(33:48)=yfc(1:16)+Ly;

        yfc(49:64)=yfc(1:16);

        yfc(65:80)=yfc(1:16)-Ly;


        %diagonals
        xfc(81:96)=xfc(1:16)+Lx;

        xfc(97:112)=xfc(1:16)-Lx;

        xfc(113:128)=xfc(1:16)-Lx;

        xfc(129:144)=xfc(1:16)+Lx;



        yfc(81:96)=yfc(1:16)+Ly;

        yfc(97:112)=yfc(1:16)+Ly;

        yfc(113:128)=yfc(1:16)-Ly;

        yfc(129:144)=yfc(1:16)-Ly;
    end
    
end

end