

%set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 16)
Volume = dX*dY*dZ;
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
for i = 1:length(thr)
    disp(i)
    M = VV*((SS^2 + thr(i)*eye(size(SS,2)))\SS*UU'); %eq(19) Thikhonov parameter needs to be determined properly.
    S = Gbig*M;
    obj(i,1) = sum(abs(S*d5-d5).^2,'all');
    obj(i,2) =  thr(i)*sum(abs(M*d5).^2,'all');
    obj(i,3) = obj(i,1)+obj(i,2);
    obj(i,4) = thr(i);
end

%%
clf
subplot(1,2,1)
loglog(thr,obj(:,1),'LineWidth',3)
[~,i1] = min(abs(thr-7e-5));
[~,i2] = min(abs(thr-8e-8));
[~,i3] = min(abs(thr-3e-10));

hold on
semilogx(thr,obj(:,3),'LineWidth',3)

plot([7e-5,7e-5+1e-11],ylim(),'k')
plot([3e-7,3e-7+1e-11],ylim(),'k')
plot([1e-9,1e-9+1e-11],ylim(),'k')
hold off
title('Tikhonov Functinal')
xlabel('\lambda')
ylabel('$J(\lambda)$')
legend('||Gm-d||_2^2','||Gm-d||_2^2+\lambda||m||_2^2','Location','southeast')
subplot(1,2,2)
loglog(thr,obj(:,2),'LineWidth',3)
hold on
plot([7e-5,7e-5+1e-11],ylim(),'k')
plot([3e-7,3e-7+1e-11],ylim(),'k')
plot([1e-9,1e-9+1e-11],ylim(),'k')

title('Regularization Term')
xlabel('\lambda')
legend('\lambda||m||_2^2')
sgtitle('Tikhonov regularisation')


%matlab2tikz('figures/tikh_RRSS.tex','width','\width','height','\height')

%% Plot solutions

   
subplot(1,3,1)
    mTilda = reshape(compute_inverse(7e-5,0,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=7e-5')
    %legend('')
subplot(1,3,2)
    mTilda = reshape(compute_inverse(3e-7,0,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet

    %caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=3e-7')
    ylabel('')
    %legend('')
subplot(1,3,3)
    mTilda = reshape(compute_inverse(1e-9,0,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=1e-9')
    ylabel('')
    %legend('')
sgtitle('Solutions with different values of k')
%matlab2tikz('figures/tikh_solutions.tex','width','\width','height','\height')

figure
    
%% Truncated SVD-regularization


thr = 1:200;
obj = zeros(length(thr),1);
for i = 1:length(thr)
    k = thr(i)
    S = Gbig*VV*eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS))*UU';
    obj(i,1) = sum(abs(S*d5-d5).^2,'all');
end

%%
clf
semilogy(thr,obj(:,1),'LineWidth',3)

title('Truncated Pseudoinverse')
xlabel('Rank of G^+ (k)')
ylabel('Residual Sum of Squares (RSS)')
legend('||Gm-d||_2^2')
[val,idx] = min(obj)
sprintf('Minimum %.3e at %d',val,thr(idx))

%matlab2tikz('figures/tsvd_RRSS.tex','width','\width','height','\height')
%% Solutions of the Truncated SVD

  %% the case when reg.parameter lambda = 0
clf
subplot(2,3,1)
    mTilda = reshape(compute_inverse(0,1,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    colorbar
    set(gca,'xtick',[])
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=1, \lambda=0')
    %legend('\lambda=8e-8, k=0')
subplot(2,3,2)
    mTilda = reshape(compute_inverse(0,5,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    %caxis([0 .25])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=5, \lambda=0')
    ylabel('')
    xlabel('')
    %legend('')
subplot(2,3,3)
    mTilda = reshape(compute_inverse(0,11,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    

    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=11, \lambda=0')
    ylabel('')
    xlabel('')
    %legend('')
    
subplot(2,3,4)
    mTilda = reshape(compute_inverse(0,20,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=23')
    %legend('')
subplot(2,3,5)
    mTilda = reshape(compute_inverse(0,41,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet

    %caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=41')
    ylabel('')
    %legend('')
subplot(2,3,6)
    mTilda = reshape(compute_inverse(0,118,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('k=118')
    ylabel('')
    %legend('')
sgtitle('Solutions for different values of k')
    
%matlab2tikz('figures/tvsd_solutions_k.tex','width','\width','height','\height')

figure
    
%% Combined Regularization 


k = 1;
Si_1 = eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS));
Mo_1 = VV*Si_1*UU';

k = 11;
Si_2 = eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS));
Mo_2 = VV*Si_2*UU';
k = 41;
Si_3 = eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS));
Mo_3 = VV*Si_3*UU';

k = 118;
Si_4 = eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS));
Mo_4 = VV*Si_4*UU';

RSS_li = [sum(abs(Gbig*Mo_1*d5-d5).^2,'all'),
          sum(abs(Gbig*Mo_2*d5-d5).^2,'all'),
          sum(abs(Gbig*Mo_3*d5-d5).^2,'all'),
          sum(abs(Gbig*Mo_4*d5-d5).^2,'all')];

thr = logspace(-14,.3,30);
%thr = [0,thr];
obj = zeros(length(thr),1);
for i = 1:length(thr)
    disp(i)
    M_1 = VV*((SS^2 + thr(i)*eye(size(SS,2)))\(SS+thr(i)*Si_1)*UU');
    M_2 = VV*((SS^2 + thr(i)*eye(size(SS,2)))\(SS+thr(i)*Si_2)*UU');
    M_3 = VV*((SS^2 + thr(i)*eye(size(SS,2)))\(SS+thr(i)*Si_3)*UU');
    M_4 = VV*((SS^2 + thr(i)*eye(size(SS,2)))\(SS+thr(i)*Si_4)*UU');
    M_4n = VV*((SS^2 + thr(i)*eye(size(SS,2)))\SS*UU');

    obj(i,1) = sum(abs(Gbig*M_1*d5-d5).^2,'all') + thr(i)*sum(abs((M_1-Mo_1)*d5).^2,'all');
    obj(i,2) = sum(abs(Gbig*M_2*d5-d5).^2,'all') + thr(i)*sum(abs((M_2-Mo_2)*d5).^2,'all');
    obj(i,3) = sum(abs(Gbig*M_3*d5-d5).^2,'all') + thr(i)*sum(abs((M_3-Mo_3)*d5).^2,'all');
    obj(i,4) = sum(abs(Gbig*M_4n*d5-d5).^2,'all') + thr(i)*sum(abs(M_4n*d5).^2,'all');
    obj(i,5) = sum(abs(Gbig*M_4*d5-d5).^2,'all') + thr(i)*sum(abs((M_4-Mo_4)*d5).^2,'all');
end

%%
clf

loglog(thr,obj(:,4),'LineWidth',3,'color',	'#D95319')
hold on

title('Tikhonov Regularization with TSVD initial value')
xlabel('\lambda')
ylabel('$J(\lambda,k)$')

loglog(thr,obj(:,1),'LineWidth',3,'color',	'#0072BD')
plot([thr(1),thr(end)], [RSS_li(1),RSS_li(1)],'--','color',	'#0072BD')
loglog(thr,obj(:,2),'LineWidth',3,'color','#7E2F8E'	)
plot([thr(1),thr(end)], [RSS_li(2),RSS_li(2)],'--','color','#7E2F8E'	)
loglog(thr,obj(:,3),'LineWidth',3,'color',	'#EDB120')
plot([thr(1),thr(end)], [RSS_li(3),RSS_li(3)],'--','color',	'#EDB120')

loglog(thr,obj(:,5),'LineWidth',3,'color',	'#4DBEEE')
plot([thr(1),thr(end)], [RSS_li(4),RSS_li(4)],'--','color',	'#4DBEEE')


plot([8e-8,8e-8+1e-11],ylim(),'k')
plot([1e-6,1e-6+1e-11],ylim(),'k')


legend('||Gm-d||_2^2+\lambda||m||_2^2','k=1','k=1 \lambda=\infty','k=11','k=11 \lambda=\infty','k=41','k=41 \lambda=\infty','k=118','k=118 \lambda=\infty','Location','southeast')

%matlab2tikz('figures/tsvd_tikh_RRSS.tex','width','\width','height','\height')
%% Solutions for combined regularisations


subplot(2,3,1)
    mTilda = reshape(compute_inverse(8e-8,1,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    caxis([0 .25])
    colorbar
    set(gca,'xtick',[])
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=8e-8, k=1')
    %legend('\lambda=8e-8, k=0')
subplot(2,3,2)
    mTilda = reshape(compute_inverse(8e-8,11,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 .25])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=8e-8, k=11')
    ylabel('')
    xlabel('')
    %legend('')
subplot(2,3,3)
    mTilda = reshape(compute_inverse(8e-8,41,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    

    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=8e-8, k=41')
    ylabel('')
    xlabel('')
    %legend('')
    
subplot(2,3,4)
    mTilda = reshape(compute_inverse(1e-6,1,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=1e-6, k=1')
    %legend('')
subplot(2,3,5)
    mTilda = reshape(compute_inverse(1e-6,11,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet

    caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=1e-6, k=11')
    ylabel('')
    %legend('')
subplot(2,3,6)
    mTilda = reshape(compute_inverse(1e-6,41,UU,SS,VV)*d5,size(Geo));
    pcolor(X,Y,abs(mTilda(:,:,7)))
    shading interp
    colormap jet
    %caxis([0 .25])
    set(gca,'ytick',[])
    colorbar
    hold on
    h = ezplot(@(x,y) (x-Xcenter).^2+(y-Ycenter).^2 - 0.015^2,[min(X), max(X), min(Y), max(Y)]);
    set(h,'color','k','LineStyle','--','LineWidth',1.5)
    hold off
    title('\lambda=1e-6, k=41')
    ylabel('')
    %legend('')
sgtitle('Solutions with different values of k')
%matlab2tikz('figures/tikh_solutions_k.tex','width','\width','height','\height')

    
