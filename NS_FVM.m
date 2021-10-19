clear,clc
close all
format compact
format long

M = 128;
rho = 10.0;
mu = 0.01;
MaxI = 900;
alpha_u = 0.050;
alpha_v = 0.050;
alpha_p = 0.005;

IW = 10;
XL = 0.5;
YL = 0.5;
NL = 100;

Vi = 1.0;
Lx = 1.0;
Ly = 1.0;

N = M;

NxU = M+1;
NyU = N+2;
NxV = M+2;
NyV = N+1;
NxP = M+2;
NyP = N+2;

dx = Lx/M;
dy = Ly/N;

hx = 1.0/dx;
hy = 1.0/dy;
vx = hx*hx;
vy = hy*hy;

data.dx = dx;
data.dy = dy;
data.hx = hx;
data.hy = hy;
data.vx = vx;
data.vy = vy;
data.U_lid = Vi;
data.rho = rho;
data.mu = mu;

u = zeros(NxU,NyU);
v = zeros(NxV,NyV);
p = zeros(NxP,NyP);

xp(1) = 0.0; xp(2:NxP-1) = 0.5*dx:dx:Lx; xp(NxP) = Lx;
yp(1) = 0.0; yp(2:NyP-1) = 0.5*dy:dy:Ly; yp(NyP) = Ly;

i = 1;
UError = 1.0;
VError = 1.0;
PError = 1.0;
SaveErrors(i,1) = i;
SaveErrors(i,2) = UError;
SaveErrors(i,3) = VError;
SaveErrors(i,4) = PError;
SaveErrors(i,5) = data.U_lid*Ly;

while ( i <= MaxI )

    u1 = u;
    v1 = v;
    p1 = p;
    
    [ us , apu ] = U_Slover( u,v,p,NxU,NyU,data,alpha_u );
    
    [ vs , apv ] = V_Slover( u,v,p,NxV,NyV,data,alpha_v );
    
    [ ps ] = P_Slover( us,apu,vs,apv,p,NxP,NyP,data,alpha_p );
    
    [ u,v,p,src ] = Correction( us,apu,vs,apv,p,ps,NxU,NyU,NxV,NyV,NxP,NyP,alpha_u,alpha_v,alpha_p,data );
    
    u2 = u;
    v2 = v;
    p2 = p;
    
    UError = (sum(sum((u2-u1).^2.0)))/(NxU*NyU);
    VError = (sum(sum((v2-v1).^2.0)))/(NxV*NyV);
    PError = (sum(sum((p2-p1).^2.0)))/(NxP*NyP);
    
    i = i+1;
    
    SaveErrors(i,1) = i;
    SaveErrors(i,2) = UError;
    SaveErrors(i,3) = VError;
    SaveErrors(i,4) = PError;
    SaveErrors(i,5) = src;
    
    if ( rem(i,IW) == 0.0 )
        Writer( u,v,p,xp,yp,data,Lx,Ly,XL,YL,NL,M,N,i );
        fprintf('Time: %4d , Xvelocity Residual: %15.12f \n',i,UError);
    end
    
end

[ s ] = Writer( u,v,p,xp,yp,data,Lx,Ly,XL,YL,NL,M,N,i );

figure(1);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,5),'b','linewidth',2);
title('Convergency of Mass','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Mass)','FontSize',14,'FontWeight','bold');

figure(2);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,2),'r','linewidth',2);
title('Convergency of Xvelocity','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Norm2(Error))','FontSize',14,'FontWeight','bold');

figure(3);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,3),'g','linewidth',2);
title('Convergency of Yvelocity','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Norm2(Error))','FontSize',14,'FontWeight','bold');

figure(4);hold on;grid on;
contourf(u',25);
title('Xvelocity','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');

figure(5);hold on;grid on;
contourf(v',25);
title('Yvelocity','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');

figure(6);hold on;grid on;
contourf(p',25);
title('Pressure','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');

save(['T',num2str(MaxI)]);








