 % Co-ordinates in nodes in global co-ordinates 
nodeheave = [   20 -18;20 -14;20 -10;20 -6;20 -4;20 -2;20 -1;20 0;   
                16 -18;16 -14;16 -10;16 -6;16 -4;16 -2;16 -1;16 0;
                12 -18;12 -14;12 -10;12 -6;12 -4;12 -2;12 -1;12 0;
                8  -18;8  -14;8  -10;8  -6;8  -4;8  -2;8  -1;8  0;
                6  -18;6  -14;6  -10;6  -6;6  -4;6  -2;6  -1;6  0;
                4  -18;4  -14;4  -10;4  -6;4  -4;4  -2;4  -1;4  0;
                3  -18;3  -14;3  -10;3  -6;3  -4;3  -2;
                2  -18;2  -14;2  -10;2  -6;2  -4;2  -2;
                1  -18;1  -14;1  -10;1  -6;1  -4;1  -2;
                0  -18;0  -14;0  -10;0  -6;0  -4;0  -2 ];

%Coordinates of nodes of individual elements in global coordinates 
eleheave =  [   1  2  10  9;2  3  11 10;3  4  12 11;4  5  13 12;5  6  14 13;6  7  15 14;7  8  16 15; %
                9  10 18 17;10 11 19 18;11 12 20 19;12 13 21 20;13 14 22 21;14 15 23 22;15 16 24 23;
                17 18 26 25;18 19 27 26;19 20 28 27;20 21 29 28;21 22 30 29;22 23 31 30;23 24 32 31;
                25 26 34 33;26 27 35 34;27 28 36 35;28 29 37 36;29 30 38 37;30 31 39 38;31 32 40 39;
                33 34 42 41;34 35 43 42;35 36 44 43;36 37 45 44;37 38 46 45;38 39 47 46;39 40 48 47;
                41 42 50 49;42 43 51 50;43 44 52 51;44 45 53 52;45 46 54 53;
                49 50 56 55;50 51 57 56;51 52 58 57;52 53 59 58;53 54 60 59;
                55 56 62 61;56 57 63 62;57 58 64 63;58 59 65 64;59 60 66 65;
                61 62 68 67;62 63 69 68;63 64 70 69;64 65 71 70;65 66 72 71  ];

totalNodes =  size(nodeheave,1);         % Total No. of nodes
ndof = totalNodes;                       % degree of freedom
numEle = size(eleheave, 1);              % No. of element
K = zeros(ndof, ndof);                   % Initializing Global Stiffness Matrix
ForceMatrix = zeros(ndof,1);             % Initializing Force Matrix


for i = 1:numEle
% Indices    
ndnumb = eleheave(i,:);
ndnumbx = nodeheave(ndnumb(1,:),1);
ndnumbz = nodeheave(ndnumb(1,:),2);

% Shape function matrix
syms eta nu 
N1 = (0.25)*(1-eta)*(1-nu);
N2 = (0.25)*(1+eta)*(1-nu);
N3 = (0.25)*(1+eta)*(1+nu);
N4 = (0.25)*(1-eta)*(1+nu);
N = [N1 N2 N3 N4]';

% Global Co-ordinates for x & z
G = [N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4]*[ndnumbx(1); ndnumbz(1); ndnumbx(2); ndnumbz(2); ndnumbx(3); ndnumbz(3); ndnumbx(4); ndnumbz(4)];

% derivative of shape functions
dXdeta = sum(diff(N.*ndnumbx,eta));
dXdnu = sum(diff(N.*ndnumbx,nu));
dZdeta = sum(diff(N.*ndnumbz,eta));
dZdnu = sum(diff(N.*ndnumbz,nu));

% Jacobian Matrix
J = [dXdeta dZdeta; dXdnu dZdnu];

% B- Matrix
B = J\[diff(N,eta) diff(N,nu)]';

ke = B'*B;
ke = int(int(ke,-1,1),-1,1);
K(ndnumb,ndnumb) = K(ndnumb,ndnumb) + ke;
end

% Free surface Boundary condition

FSN=[8  16;                         
     16 24;
     24 32;
     32 40;
     40 48];
 Kf = zeros(ndof, ndof);                   % Initializing Free surface Global Stiffness Matrix

nf = size(FSN,1);
for R=1:nf
    syms eta
    L1 = abs(nodeheave(FSN(R,1),1) - nodeheave(FSN(R,2),1));
    N1 = 1-(eta);
    N2 = eta;
    N = [N1 N2];
    FS = N'*N*(L1/2);
    FS = int(2*FS,-1,1);
    for w=1:2
        for v=1:2
            K(FSN(R,w),FSN(R,v))= K(FSN(R,w),FSN(R,v))+ FS(w,v);
        end
    end
end

%Incorporating Radiation Boundary Condition in K matrix
%Elements which comprise the Radiation Boundary
RdN =[1 2;                   
      2 3;
      3 4;
      4 5;
      5 6;
      6 7;
      7 8];

nr = size(RdN,1);
for R=1:nr
    syms nu
    L1 = abs(nodeheave(RdN(R,1),2) - nodeheave(RdN(R,2),2));
    N1 = 1-(nu);
    N2 = nu;
    N = [N1 N2];
    RN = N'*N*(L1/2);
    RdK = (-2i)*int(RN,-1,1);
    RdK = vpa(RdK,5);
    for w=1:2
        for v=1:2
            K(RdN(R,w),RdN(R,v))= K(RdN(R,w),RdN(R,v))+ RdK(w,v);
        end
    end
end

%Incorporating Body Boundary Condition into K matrix
%Elements comprising horizontal body boundary
HBBN = [46 54;               
        54 60;
        60 66;
        66 72];

% Using only horizontal body boundary 
nb1 = size(HBBN,1);
for R=1:nb1
    syms eta 
    L1 = abs(nodeheave(HBBN(R,1),1) - nodeheave(HBBN(R,2),1));
    N1 = 1-(eta);
    N2 = eta;
    N = [N1 N2];
    ny1 =  cosd(90);
    ny2 =  cosd(90);
    ny = [ny1 ny2];
    fe = (2*(L1/2)*N).*ny;
    fe = int(fe, -1,1);
    
    for w=1:2
        ForceMatrix(HBBN(R,w)) = ForceMatrix(HBBN(R,w)) + fe(w);  
    end
end

%Elements comprising vertical body boundary
VBBN =[46 47;                 
       47 48];
% Using only Vertical Body Boundary Condition
nb2 = size(VBBN,1);
for R=1:nb2
    syms nu
    L1 = abs(nodeheave(VBBN(R,1),2) - nodeheave(VBBN(R,2),2));
    N1 = 1-(nu);
    N2 = nu;
    N = [N1 N2];
    ny1 =  cosd(0);
    ny2 =  cosd(0);
    ny = [ny1 ny2];
    fe = 2*(L1/2)*N.*ny;
    fe = int(fe,-1,1);
    for w=1:2
        ForceMatrix(VBBN(R,w)) = ForceMatrix(VBBN(R,w)) + fe(w);  
    end
    
end
disp("Load Vector");
disp(ForceMatrix)


% Radiation Solution Vector Calculation
phij = vpa(K\ForceMatrix,4);
disp("Radiation Solution Matrix")
disp(phij);

density = 1030;
k = 1;                %Wave number
omega = 3.13;         %Wave frequency in Hz
g = 9.81;             %Gravitational constant

% Pressure Matrix Calculation
P = vpa(-density*i*omega*phij,4);
disp("Pressure Matrix")
disp(P);

% Horizontal Force Calculation
HorizontalForce = 0;
for i=1:nb2
    L = nodeheave(VBBN(i,1),2) - nodeheave(VBBN(i,2),2);
    n = cos(0);
    HorizontalForce = HorizontalForce + int((P(VBBN(i,1)) - P(VBBN(i,2)))*n,0,L);
end
HorizontalForce = vpa(2*HorizontalForce,4);
disp("Horizontal Force on Body Boundary")
disp(HorizontalForce);

% Added Mass Calculation
AdMB = vpa(real(HorizontalForce)/(omega*omega),5);
disp("Added Mass");
disp(AdMB);

%Added Mass Coefficient Calculation
V = 16;  
AM = vpa(AdMB/(density*V),4);
disp("Added Mass Coefficient");
disp(AM);

%Damping Coefficient Calculation
D = vpa(imag(HorizontalForce)/(density*V*omega*omega),4);
disp("Damping Coefficient");
disp(D);