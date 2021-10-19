function [ u_new,v_new,p_new,src ] = Correction( u,apu,v,apv,p,ps,NxU,NyU,NxV,NyV,NxP,NyP,alpha_u,alpha_v,alpha_p,data )

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;
vx = data.vx;
vy = data.vy;
U_lid = data.U_lid;
rho = data.rho;
mu = data.mu;

u_new = zeros(NxU,NyU);
v_new = zeros(NxV,NyV);
p_new = zeros(NxP,NyP);
source = zeros(NxP,NyP);


for i = 2:NxP-1
    for j = 2:NyP-1
        p_new(i,j) = p(i,j) + alpha_p*ps(i,j);
    end
end
i = 1;
for j = 1:NyP
    p_new(i,j) = p_new(i+1,j);
end
i = NxP;
for j = 1:NyP
    p_new(i,j) = p_new(i-1,j);
end
j = 1;
for i = 1:NxP
    p_new(i,j) = p_new(i,j+1);
end
j = NyP;
for i = 1:NxP
    p_new(i,j) = p_new(i,j-1);
end

for i = 2:NxU-1
    for j = 2:NyU-1
        up = (dy/apu(i,j))*( ps(i+1,j) - ps(i,j) );
        u_new(i,j) = u(i,j) - alpha_u*up;
    end
end
i = 1;
for j = 1:NyU
    u_new(i,j) = 0.0;
end
i = NxU;
for j = 1:NyU
    u_new(i,j) = 0.0;
end
j = 1;
for i = 1:NxU
    u_new(i,j) = -u_new(i,j+1);
end
j = NyU;
for i = 1:NxU
    u_new(i,j) = 2*U_lid-u_new(i,j-1);
end

for i = 2:NxV-1
    for j = 2:NyV-1
        vp = (dx/apv(i,j))*( ps(i,j+1) - ps(i,j) );
        v_new(i,j) = v(i,j) - alpha_v*vp;
    end
end
i = 1;
for j = 1:NyV
    v_new(i,j) = -v_new(i+1,j);
end
i = NxV;
for j = 1:NyV
    v_new(i,j) = -v_new(i-1,j);
end
j = 1;
for i = 1:NxV
    v_new(i,j) = 0.0;
end
j = NyV;
for i = 1:NxV
    v_new(i,j) = 0.0;
end

for i = 2:NxP-1
    for j = 2:NyP-1
        um = dy*rho*( u(i,j) - u(i-1,j) );
        vm = dx*rho*( v(i,j) - v(i,j-1) );
        source(i,j) = um + vm;
    end
end
src = sqrt(sum(sum(source.^2)));

end

