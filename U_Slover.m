function [ u_new , apu ] = U_Slover( u,v,p,Nx,Ny,data,WSOR )

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;
vx = data.vx;
vy = data.vy;
U_lid = data.U_lid;
rho = data.rho;
mu = data.mu;

MaxG = 10.0^-6.0;
MaxI = 25;
MinI = 10;

AP = zeros(Nx,Ny);
AE = zeros(Nx,Ny);
AW = zeros(Nx,Ny);
AN = zeros(Nx,Ny);
AS = zeros(Nx,Ny);
SP = zeros(Nx,Ny);
apu = ones(Nx,Ny);
u_new = u;
u_old = u;

A = zeros(Nx*Ny,Nx*Ny);
B = zeros(Nx*Ny,1);


De = (mu*dy*hx);
Dw = (mu*dy*hx);
Dn = (mu*dx*hy);
Ds = (mu*dx*hy);

for i = 2:Nx-1
    for j = 2:Ny-1
        ip = (j-1)*Nx+i;
        
        Fe = 0.5*rho*( u(i+1,j)   + u(i,j)   )*dy;
        Fw = 0.5*rho*( u(i-1,j)   + u(i,j)   )*dy;
        Fn = 0.5*rho*( v(i+1,j)   + v(i,j)   )*dx;
        Fs = 0.5*rho*( v(i+1,j-1) + v(i,j-1) )*dx;
        
        Dp = (p(i+1,j)-p(i,j))*dy;
        
        ae = -0.5*Fe + De;
        aw =  0.5*Fw + Dw;
        an = -0.5*Fn + Dn;
        as =  0.5*Fs + Ds;
        ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
        sp = -Dp;
        
        AP(i,j) = -ap;
        AE(i,j) =  ae;
        AW(i,j) =  aw;
        AN(i,j) =  an;
        AS(i,j) =  as;
        SP(i,j) = -sp;
        
        apu(i,j) = 1;%ap;
        
        A(ip,ip)    = AP(i,j);
        A(ip,ip+1)  = AE(i,j);
        A(ip,ip-1)  = AW(i,j);
        A(ip,ip+Nx) = AN(i,j);
        A(ip,ip-Nx) = AS(i,j);
        B(ip,1)     = SP(i,j);
        
    end
end

i = 1;
for j = 2:Ny-1
    ip = (j-1)*Nx+i;
    
    Fe = 0.5*rho*( u(i+1,j)   + u(i,j)   )*dy;
    Fw = 0.0;
    Fn = 0.5*rho*( v(i+1,j)   + v(i,j)   )*dx;
    Fs = 0.5*rho*( v(i+1,j-1) + v(i,j-1) )*dx;
    
    ae = -0.5*Fe + De;
    aw =  0.0;
    an = -0.5*Fn + Dn;
    as =  0.5*Fs + Ds;
    ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
    
    apu(i,j) = 1;%ap;
    
end


i = Nx;
for j = 2:Ny-1
    ip = (j-1)*Nx+i;
    
    Fe = 0.0;
    Fw = 0.5*rho*( u(i-1,j)   + u(i,j)   )*dy;
    Fn = 0.5*rho*( v(i+1,j)   + v(i,j)   )*dx;
    Fs = 0.5*rho*( v(i+1,j-1) + v(i,j-1) )*dx;
    
    ae = -0.0;
    aw =  0.5*Fw + Dw;
    an = -0.5*Fn + Dn;
    as =  0.5*Fs + Ds;
    ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
    
    apu(i,j) = 1;%ap;
    
end

GError = 1.0;
iter = 1;
while ( GError > MaxG )
    
    i = 1;
    for j = 1:Ny
        u_new(i,j) = 0;
    end
    i = Nx;
    for j = 1:Ny
        u_new(i,j) = 0;
    end
    
    j = 1;
    for i = 1:Nx
        u_new(i,j) = -u_new(i,j+1);
    end
    j = Ny;
    for i = 1:Nx
        u_new(i,j) = 2*U_lid-u_new(i,j-1);
    end
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            Be = AE(i,j)*u_old(i+1,j);
            Bw = AW(i,j)*u_new(i-1,j);
            Bn = AN(i,j)*u_old(i,j+1);
            Bs = AS(i,j)*u_new(i,j-1);
            u_new(i,j) = (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            u_new(i,j) = WSOR*u_new(i,j)+(1.0-WSOR)*u_old(i,j);
        end
    end
    
    GError = (sum(sum((u_new-u_old).^2.0)))/(Nx*Ny);
    
    u_old = u_new;
    
    iter = iter+1;
    
    if( iter > MaxI )
        break;
    end
    if( iter < MinI )
        GError = 1.0;
    end
    
end


end

