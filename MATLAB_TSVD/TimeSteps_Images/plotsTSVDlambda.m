
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


%lambda = 0;

% this lambda corresponds to the sigma(136)
lambda = 1.7333*10^-8

 % lambda = 1;
  
 
thr = 1:136;
obj = zeros(length(thr),3);

timemoments = 5;

rhs = zeros(5,length(d1));
Residual  = zeros(5, length(thr));
Regterm  = zeros(5, length(thr));
Tikfunc  = zeros(5, length(thr));

rhs(1,:) =  d1;
rhs(2,:) =  d2;
rhs(3,:) =  d3;
rhs(4,:) =  d4;
rhs(5,:) =  d5;


for k = 1:5

	d =  rhs(k,:)';


for i = 1:length(thr)
	  j = thr(i);
  %  S = Gbig*VV*eye(size(SS,1),j)*eye(j,size(SS,2))*diag(1./diag(SS))*UU';
%here we compute matrix M using TSVD for index j and reg. parameter lambda
  M = compute_inverse(lambda,j,UU,SS,VV);
  
  S = Gbig*M;

  % here we compute residuals
    obj(i,1) = sum(abs(S*d-d).^2,'all');
    obj(i,2) =  lambda*sum(abs(M*d).^2,'all');
    obj(i,3) = obj(i,1)+obj(i,2);

end

%save values of the residual
Residual(k,:) = obj(:,1);
Regterm(k,:)  = obj(:,2);
Tikfunc(k,:) =  obj(:,3);

end

% plotting dependence of the residual from truncation parameter k in TSVD
%****************************************************************************
figure

for k = 1:5
	  subplot(5,1,k)

	%  loglog(thr,Residual(k,:),'LineWidth',3)
	 semilogy(thr,Residual(k,:),'LineWidth',3)
	  
	  xlabel('Rank of G^+ (k)')
	  ylabel('R(O_k)')
	  legend('R(O_k) = ||G_k O_k - d||_2^2')
	  title([' d_i, i=',num2str(k)])

   [val,idx] = min(Residual(k,:))
   sprintf('Minimum of residual is %.3e at  k= %d',val,thr(idx))
	  MinResid(k) = val; 	  
	  end

  
    figure

        for k = 1:5

	     semilogy(thr,Residual(k,:),'LineWidth',3)
	     hold on
	     end
		      
		  xlabel('k');

	    
legend('t=2 min','t=4 min','t=6 min','t=8 min','t=10 min');
		      
title(['|| G_k O_k - d ||_2^2, \lambda=',num2str(lambda)]);


hold off
	  
	  %*****************************************************************
	  figure

   semilogy(MinResid, "--ks", "MarkerSize", 7, "MarkerFaceColor", "m")


  %str_xlabel = ['lambda =', num2str(lambda)];

%title(['Minimal residual for lambda= ',num2str(lambda),', ',str_xlabel])
 
title([' ||G_k O_k - d||_2^2,  k=136, \lambda= ',num2str(lambda)])
  legend('Time moments t = 2,4,6,8,10 (min)')




%***********************************************************************
% plotting reg.term depending on the truncated value  k
%***********************************************************************

figure

for k = 1:5
	  subplot(5,1,k)

	 semilogy(thr,Regterm(k,:),'LineWidth',3)
	  
	  xlabel('Rank of G^+ (k)')
%	  ylabel('\lambda ||O_k||_2^2')
	  legend('\lambda ||O_k||_2^2')
	  title([' d_i, i=',num2str(k)])

   [val,idx] = min(Regterm(k,:))
   sprintf('Minimum of the reg.term is %.3e at  k= %d',val,thr(idx))
	  MinReg(k) = val; 	  
	  end

 %**************************************************************
	  figure

   semilogy(MinReg, "--ks", "MarkerSize", 7, "MarkerFaceColor", "m") 
 title([' \lambda ||O_k||_2^2 for k=136 and \lambda= ',num2str(lambda)])

 legend('Time moments t = 2,4,6,8,10 (min)')

	    
%***********************************************************************
% plotting Tikh.func = Minim. residual + reg.term depending on the truncated value  k
%***********************************************************************

figure

for k = 1:5
	  subplot(5,1,k)

	 semilogy(thr,Tikfunc(k,:),'LineWidth',3)
	  
	  xlabel('Rank of G^+ (k)')

	  legend('|| G_k O_k - d ||_2^2 + \lambda ||O_k||_2^2')
	  title([' d_i, i=',num2str(k),' \lambda=',num2str(lambda)])

   [val,idx] = min(Tikfunc(k,:))
   sprintf('Minimum of the functional is %.3e at  k= %d',val,thr(idx))
	  MinTik(k) = val; 	  
	  end

	  figure

   semilogy(MinTik, "--ks", "MarkerSize", 7, "MarkerFaceColor", "m") 
 title(['|| G_k O_k - d ||_2^2 + \lambda ||O_k||_2^2 for k=136 and \lambda= ',num2str(lambda)])

	    legend('Time moments t = 2,4,6,8,10 (min)')
	    
%*************************************************************************	    
    figure

        for k = 1:5

	     semilogy(thr,Tikfunc(k,:),'LineWidth',3)
	     hold on
	     end
		      
		  xlabel('k');

	    
legend('t=2 min','t=4 min','t=6 min','t=8 min','t=10 min');
		      
title(['|| G_k O_k - d ||_2^2 + \lambda ||O_k||_2^2, \lambda=',num2str(lambda)]);


hold off

%****************************************************************************
