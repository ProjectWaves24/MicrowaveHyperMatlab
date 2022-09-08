
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



thr = 1:256;
obj = zeros(length(thr),1);

timemoments = 5;

rhs = zeros(5,length(d1));
Residual  = zeros(5, length(thr));

rhs(1,:) =  d1;
rhs(2,:) =  d2;
rhs(3,:) =  d3;
rhs(4,:) =  d4;
rhs(5,:) =  d5;


%in this program the reg.parameter lambda = 0
for k = 1:5

	d =  rhs(k,:)';


for i = 1:length(thr)
	  j = thr(i);
    S = Gbig*VV*eye(size(SS,1),j)*eye(j,size(SS,2))*diag(1./diag(SS))*UU';
    obj(i,1) = sum(abs(S*d-d).^2,'all');
end

%save values of the residual
Residual(k,:) = obj(:,1);

end



% plotting dependence of the residual from truncation parameter k in TSVD
figure

for k = 1:5
	  subplot(5,1,k)

	%  loglog(thr,Residual(k,:),'LineWidth',3)
	 semilogy(thr,Residual(k,:),'LineWidth',3)
	  
	  xlabel('Rank of G^+ (k)')
	  ylabel('R(m)')
	  legend('R(m) = ||Gm-d||_2^2')
	  title([' d(i), i=',num2str(k)])

   [val,idx] = min(Residual(k,:))
sprintf('Minimum of residual is %.3e at  k= %d',val,thr(idx))
	  end


	  
