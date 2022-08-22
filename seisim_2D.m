% Ruhul Firdaus
% April 2022

clear; close all; clc;
load model3.mat
c = Param_2;

% reframe
dz = 5; % meter
zz = 5000; % meter
z = 0:dz:zz;
Nz = length(z);

dx = 50; % meter
xx = 15000; % meter
x = 0:dx:xx;
Nx = length(x);

line_start_at = 5000;

close all

%% Simulation param
% source depth
z_s = 0;
z_s = round(z_s/dz)+1;
z_r = 0;
z_r = round(z_r/dz)+1;

% recompute dt for stability criteria
c_ref = 1.1*max(max(c)); % c max
CFL = 1; 
dt = CFL*min(dx,dz)/c_ref; % CFL may not exceed 1

% time increment
tt = 5; %detik
t = 0:dt:tt;
Nt = length(t);

% remodel Ricker source for new dt
tt_r = 0.5; %detik
t_r = -tt_r/2:dt:tt_r/2;

f_Scale = tt_r/10; % atur freq wavelet

g = normpdf(t_r,0,f_Scale);
r = (g(3:end)-g(1:end-2))/2/dt;
r = (r(3:end)-r(1:end-2))/2/dt;
r = -[0 0 r 0 0];
r = r/max(abs(r));
figure(80); plot(t_r,r); grid on; grid minor; title 'Source : Ricker'
xlabel 't(s)'; xlim([t_r(1) t_r(end)]); ylim([-1.2 1.2]);

% Fourier transform analysis
td_s = r;
fta

% redefine t as a consequence of zero phase wavelet
t = t - tt_r/2;


%% 
% resize model parameter
c = resizem(c,[Nz Nx],'nearest');

f = figure(303);f.WindowState = 'maximized';
subplot(1,2,1);
imagesc(x,z,c); title 'velocity model'; xlabel 'x(m)'; ylabel 'z(m)';
colormap 'bone'; hold on;

%% Desain parameter akuisisi
close all;
N_channel = 50; % genap
g_int = dx; % meter
N_src = 3;
s_int = 25*g_int; % meter
near_offset = 0;
% sampling_rate = dt; % second
% record length = tt/dt; % second
% far_offset = (N_channel-1)*g_int + near_offset;
% spread_length = far_offset;
% line_spread = ;

spread_type = 'fixed_recs'; %'off_end'; 'split'; 'fixed_recs'

s_loc=[]; r_loc=[];
if strcmp(spread_type,'off_end')
    for i = 1:N_src
        s_temp = line_start_at+(i-1)*s_int;
        for j = 1:N_channel
            s_loc = [s_loc s_temp]; % meter
            g_temp = s_temp+near_offset + (j-1)*g_int;
            r_loc = [r_loc g_temp]; % meter
        end
    end
    
elseif strcmp(spread_type, 'split')
    s_loc=[]; r_loc=[];
    for i = 1:N_src
        s_temp1 = line_start_at+(i-1)*s_int;
        s_temp2 = line_start_at+(i-1)*s_int;
        for j = 1:N_channel/2
            s_loc = [s_temp2 s_loc s_temp1]; % meter
            g_temp1 = s_temp1+near_offset + (j-1)*g_int;
            g_temp2 = s_temp2-near_offset - (j-1)*g_int;
            r_loc = [g_temp2 r_loc g_temp1]; % meter
        end
    end
elseif strcmp(spread_type,'fixed_recs')
    for i = 1:N_src
        s_temp = line_start_at+(i-1)*s_int;
        for j = 1:N_channel
            s_loc = [s_loc s_temp]; % meter
            g_temp = line_start_at + near_offset + (j-1)*g_int;
            r_loc = [r_loc g_temp]; % meter
        end
    end

end

figure(653);plot((r_loc+s_loc)/2,s_loc,'og');
hold on;
plot(r_loc,s_loc,'.b'); xlabel 'x_{rec} (m)'; ylabel 'x_{src} (m)'; 
grid on; grid minor;ylim([x(1) x(end)]);xlim([x(1) x(end)]);

r_loc = reshape(r_loc,[N_channel N_src]); % meter
s_loc = reshape(s_loc,[N_channel N_src]); % meter

% converting g_loc and s_loc to grid position
r_loc = round(r_loc/dx)+1; % grid point
s_loc = round(s_loc/dx)+1; % grid point

%%
% set recording format
R = zeros(Nt,N_channel,N_src);
G = zeros(Nt,N_channel,N_src);
% U = zeros(Nz,Nx,Nt,N_src);

for is = 1:N_src
    nshot = is
    
    % initial condition
    u = zeros(Nz,Nx,Nt);
    u = single(u);

    for it = 2:Nt-1
        progress = round(it/Nt*100)
        % source position
        u(z_s,s_loc(1,is),1:length(r)) = r;

        % Reverse time
%         u(1,g_loc,:) = flipud(permute(R,[3 2 1]));
%         temp = flipud(permute(U(1,:,:),[3 2 1]));
%         u(1,:,:) = permute(temp,[3 2 1]);
        
        for iz = 2:Nz-1
            for ix = 2:Nx-1
                kx = (dt/dx*c(iz,ix));
                kz = (dt/dz*c(iz,ix));
                u(iz,ix,it+1) = 2*u(iz,ix,it) - u(iz,ix,it-1) +...
                    kx^2 * ( u(iz,ix+1,it) - 2 * u(iz,ix,it) + u(iz,ix-1,it) )+...
                    kz^2 * ( u(iz+1,ix,it) - 2 * u(iz,ix,it) + u(iz-1,ix,it) );
            end
        end

        % Absorbing Boundary Condition
        u(:,  1, it+1) = u(:,  1, it) + (c(iz, 1)*dt/dx) * ( u(:, 2, it) - u(:, 1,   it) ); %centered difference in x
        u(:, Nx, it+1) = u(:, Nx, it) - (c(iz,Nx)*dt/dx) * ( u(:, Nx,it) - u(:, Nx-1,it) ); %centered difference in x
        u(1 , :, it+1) = u(1 , :, it) + (c(1, ix)*dt/dz) * ( u(2 , :,it) - u(1   , :,it) ); %centered difference in x
        u(Nz, :, it+1) = u(Nz, :, it) - (c(Nz,ix)*dt/dz) * ( u(Nz, :,it) - u(Nz-1, :,it) ); %centered difference in x

    end
 
    for ir = 1:N_channel
        R_temp = u(z_r,r_loc(ir,is),:);
        R_temp = permute(R_temp,[3 1 2]);
        R(:,ir,is) = R_temp;
        G_temp = u(2*round(Nz/5),r_loc(ir,is),:);
        G_temp = permute(G_temp,[3 1 2]);
        G(:,ir,is) = G_temp;
    end
    
    % pre-set figure size
    f = figure(303);f.WindowState = 'maximized';

    mm = 50*mean(mean(c));
    cm = [min(min(c)) max(max(c))];
    dit = 50;
    
    show = 1;
    if show==1
        figure(303);clf
        subplot(1,2,2)
        imagesc(x,t,permute(u(3,:,:),[3 2 1]));caxis([-0.01 0.01]);colormap 'bone'
        xlabel 'x(m)'; ylabel 't(s)'; title(['Shot ',num2str(is)])

        subplot(1,2,1);

        for it = 2:dit:Nt-1
            imagesc(x,z,c+mm*u(:,:,it+1)); colormap 'bone'; hold on;
            grid on; grid minor
            plot((r_loc(:,is)-1)*dx,(z_r-1)*dz,'.b','MarkerSize',5); 
            plot((s_loc(:,is)-1)*dx,(z_s-1)*dz,'.r','MarkerSize',5);
            title 'velocity model with perturbation';
            title(['t = ',num2str(t(it+1)),' s']);
            xlabel 'x(m)'; ylabel 'z(m)'
            caxis(cm);
            pause(0.0001)
        end
    end
    
    U{is} = u;

end

RR = reshape(R,[Nt N_src*N_channel]);
GG = reshape(G,[Nt N_src*N_channel]);

%%
show = 1;
% pre-set figure size
f = figure(303);f.WindowState = 'maximized';

mm = 50*mean(mean(c));
cm = [min(min(c)) max(max(c))];
dit = 50;
dis = 1;

if show==1
    for is = 1:dis:N_src
        is
        u = U{is};
        figure(303);clf
        subplot(1,2,2)
        imagesc(x,t,permute(u(3,:,:),[3 2 1]));caxis([-0.01 0.01]);colormap 'bone'
        xlabel 'x(m)'; ylabel 't(s)'; title(['Shot ',num2str(is)])

        subplot(1,2,1);
        
        for it = 2:dit:Nt-1
            imagesc(x,z,c+mm*u(:,:,it+1)); colormap 'bone'; hold on;
            grid on; grid minor
            plot((r_loc(:,is)-1)*dx,(z_r-1)*dz,'.b','MarkerSize',5); 
            plot((s_loc(:,is)-1)*dx,(z_s-1)*dz,'.r','MarkerSize',5);
            title 'velocity model with perturbation';
            title(['t = ',num2str(t(it+1)),' s']);
            xlabel 'x(m)'; ylabel 'z(m)'
            caxis(cm);
            pause(0.0001)
        end
    end
end

%% resample
t_old = t;
dt = 1e-3;
% time increment
t = 0:dt:tt;
Nt = length(t);

% d1 = 1;
% for d1 = d1:size(u,1)
%     d1
%     for d2 = 1:size(u,2)
%         uu(d1,d2,:)=interp1(t_old,squeeze(u(d1,d2,:)),t);
%     end
% end
% u = uu; clear uu

d2 = 1;
for d2 = d2:size(R,2)
    for d3 = 1:size(R,3)
        Rr(:,d2,d3)=interp1(t_old,squeeze(R(:,d2,d3)),t);
    end
end
R = Rr; clear Rr

d2 = 1;
for d2 = d2:size(G,2)
    for d3 = 1:size(G,3)
        Gg(:,d2,d3)=interp1(t_old,squeeze(G(:,d2,d3)),t);
    end
end
G = Gg; clear Gg

d2 = 1;
for d2 = d2:size(RR,2)
    for d3 = 1:size(RR,3)
        RRr(:,d2,d3)=interp1(t_old,squeeze(RR(:,d2,d3)),t);
    end
end
RR = RRr; clear RRr

d2 = 1;
for d2 = d2:size(GG,2)
    for d3 = 1:size(GG,3)
        GGg(:,d2,d3)=interp1(t_old,squeeze(GG(:,d2,d3)),t);
    end
end
GG = GGg; clear GGg

%%
figure(403);
imagesc(x,t,R(:,:,1));caxis([-0.01 0.01]);colormap 'bone'
xlabel 'x(m)'; ylabel 't(s)'

figure(503);
imagesc(x,t,RR);caxis([-0.01 0.01]);colormap 'bone'
xlabel 'x(m)'; ylabel 't(s)'

figure(404);
imagesc(x,t,G(:,:,1));caxis([-0.01 0.01]);colormap 'bone'
xlabel 'x(m)'; ylabel 't(s)'

figure(504);
imagesc(x,t,GG);caxis([-0.01 0.01]);colormap 'bone'
xlabel 'x(m)'; ylabel 't(s)'
