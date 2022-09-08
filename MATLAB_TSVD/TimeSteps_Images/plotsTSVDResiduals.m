
close all

%set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 16)
%Volume = dX*dY*dZ;
%%
% Requires 
% Precomputed_Inverse_Scattering_Solution_M.m
%   - UU, SS, VV = SVD(G)
% ImageWrapper_TimeSeries.m
%   - DeltaS
%   - NumAnt
t = 4;
d1 = reshape(DeltaS(:,:,1),NumAnt^2,1);
d2 = reshape(DeltaS(:,:,2),NumAnt^2,1);
d3 = reshape(DeltaS(:,:,3),NumAnt^2,1);
d4 = reshape(DeltaS(:,:,4),NumAnt^2,1);
d5 = reshape(DeltaS(:,:,5),NumAnt^2,1);

%%

thr = logspace(-14,3,30);
%thr = [0,thr];
obj = zeros(length(thr),1);

timemoments = 5;

rhs = zeros(5,length(d1));
Tikhfunc  = zeros(5, length(thr));

rhs(1,:) =  d1;
rhs(2,:) =  d2;
rhs(3,:) =  d3;
rhs(4,:) =  d4;
rhs(5,:) =  d5;



figure

for k = 1:5

	d =  rhs(k,:)';


for i = 1:length(thr)
 %   disp(i)
	  truncation = 9;
 % here we compute matrix using TSVD
 M=  compute_inverse(thr(i),truncation,UU,SS,VV);
	  
  %  M = VV*((SS^2 + thr(i)*eye(size(SS,2)))\SS*UU');
    
    S = Gbig*M;
 
% Compute L_2 norm for residual  R(m) = || Gm -d||_2^2
  obj(i,1) = sum(abs(S*d - d).^2,'all');
 

    obj(i,2) =  thr(i)*sum(abs(M*d).^2,'all');
 
   % obj(i,3) = obj(i,1)+obj(i,2);
    obj(i,3) = obj(i,1);
    obj(i,4) = thr(i);
	       
end

%save values of Tikh.func.
Tikhfunc(k,:) = obj(:,1);

%save values of reg.term
Regterm(k,:) =  obj(:,2);

% plotting dependence of the residual from reg.parameter

subplot(5,1,k)
    loglog(thr,Tikhfunc(k,:),'LineWidth',3)
xlabel('\lambda')
 ylabel('R(O)')
 legend('R(O) = ||GO_9 - d||_2^2')
 title([' d(i), i=',num2str(k)])

 
end
           


%%

%plotting dependence of the reg.term from the reg.parameter

  figure

for k = 1:5  
	  subplot(5,1,k)
	  loglog(thr,Regterm(k,:),'LineWidth',3)
	  xlabel('\lambda')

	  legend('\lambda || O_9 ||_2^2')
	  title(['d(i), i=',num2str(k)])
end
 
