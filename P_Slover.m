function [ p_new ] = P_Slover( us,apu,vs,apv,p,Nx,Ny,data,WSOR )

SType = 2;

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;
vx = data.vx;
vy = data.vy;
U_lid = data.U_lid;
rho = data.rho;
mu = data.mu;

MaxG = 10.0^-8.0;
MaxI = 7500;
MinI = 6000;

AP = zeros(Nx,Ny);
AE = zeros(Nx,Ny);
AW = zeros(Nx,Ny);
AN = zeros(Nx,Ny);
AS = zeros(Nx,Ny);
SP = zeros(Nx,Ny);

p_new = p;
p_old = p;

A = zeros(Nx*Ny,Nx*Ny);
B = zeros(Nx*Ny,1);

for i = 2:Nx-1
    for j = 2:Ny-1
        ip = (j-1)*Nx+i;
        
        um = dy*rho*( us(i,j) - us(i-1,j) );
        vm = dx*rho*( vs(i,j) - vs(i,j-1) );
        
        ae = rho*(1.0/vy)/apu(i,j);
        aw = rho*(1.0/vy)/apu(i-1,j);
        an = rho*(1.0/vx)/apv(i,j);
        as = rho*(1.0/vx)/apv(i,j-1);
        ap = ae + aw + an + as;
        
        sp = - um - vm;
        
        AP(i,j) = -ap;
        AE(i,j) =  ae;
        AW(i,j) =  aw;
        AN(i,j) =  an;
        AS(i,j) =  as;
        SP(i,j) = -sp;
        
        A(ip,ip)    = AP(i,j);
        A(ip,ip+1)  = AE(i,j);
        A(ip,ip-1)  = AW(i,j);
        A(ip,ip+Nx) = AN(i,j);
        A(ip,ip-Nx) = AS(i,j);
        B(ip,1)     = SP(i,j);
        
    end
end

AP(2,2) = 1.0;
A(Nx+2,Nx+2) = AP(2,2);
B(Nx+2,1) = 0.0;

if ( SType == 1 )
    
    i = 1;
    for j = 1:Ny
        ip = (j-1)*Nx+i;
        A(ip,ip) = 1.0;
        A(ip,ip+1) = -1.0;
        B(ip,1) = 0.0;
    end
    i = Nx;
    for j = 1:Ny
        ip = (j-1)*Nx+i;
        A(ip,ip) = 1.0;
        A(ip,ip-1) = -1.0;
        B(ip,1) = 0.0;
    end
    
    j = 1;
    for i = 2:Nx-1
        ip = (j-1)*Nx+i;
        A(ip,ip) = 1.0;
        A(ip,ip+Nx) = -1.0;
        B(ip,1) = 0.0;
    end
    j = Ny;
    for i = 2:Nx-1
        ip = (j-1)*Nx+i;
        A(ip,ip) = 1.0;
        A(ip,ip-Nx) = -1.0;
        B(ip,1) = 0.0;
    end
    
    Xp = A\B;
    
    for i = 1:Nx
        for j = 1:Ny
            ip = (j-1)*Nx+i;
            p_new(i,j) = Xp(ip);
        end
    end
    
else
    
    GError = 1.0;
    iter = 1;
    while ( GError > MaxG )
        
        i = 1;
        for j = 1:Ny
            p_new(i,j) = p_new(i+1,j);
        end
        i = Nx;
        for j = 1:Ny
            p_new(i,j) = p_new(i-1,j);
        end
        
        j = 1;
        for i = 1:Nx
            p_new(i,j) = p_new(i,j+1);
        end
        j = Ny;
        for i = 1:Nx
            p_new(i,j) = p_new(i,j-1);
        end
        
        for i = 2:Nx-1
            for j = 2:Ny-1
                Be = AE(i,j)*p_old(i+1,j);
                Bw = AW(i,j)*p_new(i-1,j);
                Bn = AN(i,j)*p_old(i,j+1);
                Bs = AS(i,j)*p_new(i,j-1);
                p_new(i,j) = (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            end
        end
        
        GError = (sum(sum((p_new-p_old).^2.0)))/(Nx*Ny);
        
        p_old = WSOR*p_new+(1.0-WSOR)*p_old;
        
        iter = iter+1;
        
        if( iter > MaxI )
            break;
        end
        if( iter < MinI )
            GError = 1.0;
        end
        
    end
    
end

end

