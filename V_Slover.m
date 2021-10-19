function [ v_new , apv ] = V_Slover( u,v,p,Nx,Ny,data,WSOR )

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
apv = ones(Nx,Ny);
v_new = v;
v_old = v;

A = zeros(Nx*Ny,Nx*Ny);
B = zeros(Nx*Ny,1);

De = (mu*dy*hx);
Dw = (mu*dy*hx);
Dn = (mu*dx*hy);
Ds = (mu*dx*hy);

for i = 2:Nx-1
    for j = 2:Ny-1
        ip = (j-1)*Nx+i;
        
        Fe = 0.5*rho*( u(i,j+1)   + u(i,j)   )*dy;
        Fw = 0.5*rho*( u(i-1,j+1) + u(i-1,j) )*dy;
        Fn = 0.5*rho*( v(i,j+1)   + v(i,j)   )*dx;
        Fs = 0.5*rho*( v(i,j-1)   + v(i,j)   )*dx;
        
        Dp = (p(i,j+1)-p(i,j))*dx;
        
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
        
        apv(i,j) = 1;%ap;
        
        A(ip,ip)    = AP(i,j);
        A(ip,ip+1)  = AE(i,j);
        A(ip,ip-1)  = AW(i,j);
        A(ip,ip+Nx) = AN(i,j);
        A(ip,ip-Nx) = AS(i,j);
        B(ip,1)     = SP(i,j);
        
    end
end

j = 1;
for i = 2:Nx-1
    ip = (j-1)*Nx+i;
    
    Fe = 0.5*rho*( u(i,j+1)   + u(i,j)   )*dy;
    Fw = 0.5*rho*( u(i-1,j+1) + u(i-1,j) )*dy;
    Fn = 0.5*rho*( v(i,j+1)   + v(i,j)   )*dx;
    Fs = 0.0;
    
    ae = -0.5*Fe + De;
    aw =  0.5*Fw + Dw;
    an = -0.5*Fn + Dn;
    as = 0.0;
    ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
    
    apv(i,j) = 1;%ap;
    
end

j = Ny;
for i = 2:Nx-1
    ip = (j-1)*Nx+i;
    
    Fe = 0.5*rho*( u(i,j+1)   + u(i,j)   )*dy;
    Fw = 0.5*rho*( u(i-1,j+1) + u(i-1,j) )*dy;
    Fn = 0.0;
    Fs = 0.5*rho*( v(i,j-1)   + v(i,j)   )*dx;
    
    ae = -0.5*Fe + De;
    aw =  0.5*Fw + Dw;
    an = 0.0;
    as =  0.5*Fs + Ds;
    ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
    
    apv(i,j) = 1;%ap;
    
end

GError = 1.0;
iter = 1;
while ( GError > MaxG )
    
    i = 1;
    for j = 1:Ny
        v_new(i,j) = -v_new(i+1,j);
    end
    i = Nx;
    for j = 1:Ny
        v_new(i,j) = -v_new(i-1,j);
    end
    
    j = 1;
    for i = 1:Nx
        v_new(i,j) = 0.0;
    end
    j = Ny;
    for i = 1:Nx
        v_new(i,j) = 0.0;
    end
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            Be = AE(i,j)*v_old(i+1,j);
            Bw = AW(i,j)*v_new(i-1,j);
            Bn = AN(i,j)*v_old(i,j+1);
            Bs = AS(i,j)*v_new(i,j-1);
            v_new(i,j) = (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            v_new(i,j) = WSOR*v_new(i,j)+(1.0-WSOR)*v_old(i,j);
        end
    end
    
    GError = (sum(sum((v_new-v_old).^2.0)))/(Nx*Ny);
    
    v_old = v_new;
    
    iter = iter+1;
    
    if( iter > MaxI )
        break;
    end
    if( iter < MinI )
        GError = 1.0;
    end
    
end


end

