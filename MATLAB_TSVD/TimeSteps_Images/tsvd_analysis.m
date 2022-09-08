close all
[~,dirname] = fileparts(pwd);


% run first the program Precomputed_Inverse_Scattering_SOlution_M.m

% x-axis to plot singular varlues

x = linspace(0, 255, 256);


sdiag = diag(SS);

figure

% plotting all singular values
plot(x, diag(SS), 'LineWidth',2);
xlabel('k')


legend('singular values  of \Sigma')
 saveas(gcf,append(dirname,sprintf('Sigmaall.png')))

%**************************************************************

figure

% plotting singular values
plot(x(1:136), sdiag(1:136), 'LineWidth',2);
xlabel('k')

legend('first 136 singular values of \Sigma')

 saveas(gcf,append(dirname,sprintf('Sigma136.png')))
%************************************************************

figure


% plotting  1/ Sigma
plot(x, 1./diag(SS),'LineWidth',2);
xlabel('k')

legend('1./\Sigma')

%title(' \lambda=0')

saveas(gcf,append(dirname,sprintf('InvSigma.png')))


figure
 
%******************************************************
% plotting of 1/Sigma for lambda=0
  % this test explains why for small trunkated values we get good reconstruction
  % Accordingly Corollary for Theorem 5 we get error 1/sigma_r *norm |d|
  % Here we see that with increasing of trunkation number k increases also the error.
plot(x(1:136), 1./sdiag(1:136),'LineWidth',2);

xlabel('k')

legend('first 136 values  of 1./\Sigma')
%title(' \lambda=0')
 saveas(gcf,append(dirname,sprintf('InvSigma136.png')))
%************************************************************
figure

lambda = sdiag(1);

% plotting singular values
plot(x(1:136), sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda),'r --', 'LineWidth',2);

hold on


lambda = sdiag(136);

plot(x(1:136), sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda),'b -','LineWidth',2);
xlabel('k')


legend('\sigma_{k+1}/(\sigma_{k+1}^2 + \sigma_1)','\sigma_{k+1}/(\sigma_{k+1}^2 +  \sigma_{136})')

title('first 136 values of \sigma_{k+1}/(\sigma_{k+1}^2 + \lambda)')

 saveas(gcf,append(dirname,sprintf('func136.png')))
 
 
sprintf('Maximal indexes for Sigma > sqrt(lambda), sigma(136)=lambda')
 find(sdiag > sqrt(lambda))' 

figure

lambda = sdiag(1);

plot(x(1:136), sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda),'LineWidth',2);
title('first 136 values of \sigma_{k+1}/(\sigma_{k+1}^2 + \lambda)')
legend('\sigma_{k+1}/(\sigma_{k+1}^2 + \sigma_1)','\sigma_{k+1}/(\sigma_{k+1}^2 +  \sigma_{1})')

saveas(gcf,append(dirname,sprintf('func1.png')))
 
 
%****************************************************************************
figure

lambda = sdiag(136);

% plotting singular values
plot(x(1:136), sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda),'LineWidth',2);
xlabel('k')


legend('\sigma_{k+1}/(\sigma_{k+1}^2 + \sigma_{136})')

title('first 136 values of \sigma_{k+1}/(\sigma_{k+1}^2 + \lambda)')

%************************************************************************

figure


  
% plotting singular values for small values of the reg.parameter


%discr_lambda = linspace(1,136,6);
discr_lambda = linspace(1,256,6)


  
for i=1:length(discr_lambda)
    
	lambda = sdiag(discr_lambda(i));

for j=1:136
 % rel_error_optimalr(j) = sdiag(136)/(sdiag(136)*sdiag(136) + lambda);
  rel_error(i,j) = sdiag(j)/(sdiag(j)*sdiag(j) + lambda);
end

sprintf('Maximal indexes for Sigma > sqrt(lambda)')
 find(sdiag > sqrt(lambda))' 
 
% sprintf('Indexes for Sigma  when Sigma_k = lambda')
%find(sdiag == lambda)' 

  % here is computation of the bound in the Theorem for large lambda > sigma_i
  
%rel_error_optimalr(i,:) = sdiag(136)./(sdiag(136)*sdiag(136) + lambda);
%rel_error(i,:) = sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda);



subplot(6,1,i)
plot(x(1:136), rel_error(i,:),'LineWidth',2);

%hold on
%plot(x(1:136),rel_error_optimalr,'r --', 'LineWidth',2);

  title(['\lambda=',num2str(lambda)])
 
 %  hold off
  
end


xlabel('k')

 % title(['d(i), i=',num2str(k)])

%title('Optimal \lambda=\sigma_{k+1}')

  
 saveas(gcf,append(dirname,sprintf('optimallambda.png')))

  
%****************  plotting  similar figures for large values of lambda
%*****************************************************************************


figure


  
% plotting singular values
%discr_lambda = linspace(1,136,4);
%discr_lambda = linspace(1,256,6)

    % here we produce large values of lambda from 10^1 to 10^6
discr_lambda = logspace(1,6,6);


for i=1:length(discr_lambda)
	lambda = discr_lambda(i);


  % here is computation of the bound in the Theorem for large lambda > sigma_i
  
for j=1:136
  rel_error_optimalr(j) = sdiag(136)/(sdiag(136)*sdiag(136) + lambda);
  rel_error(i,j) = sdiag(j)/(sdiag(j)*sdiag(j) + lambda);
end

%rel_error(i,:) = sdiag(1:136)./(sdiag(1:136).*sdiag(1:136) + lambda);


subplot(6,1,i)
plot(x(1:136), rel_error(i,:),'LineWidth',2);


  title(['\lambda=',num2str(lambda)])
  
end


xlabel('k')
