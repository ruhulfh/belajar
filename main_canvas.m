clear;close all; clc;
Fg = figure(7485);
x_min = 0;
x_max = 5000;
z_min = 0;
z_max = 1000;
xlim([x_min x_max])
ylim([z_min z_max])
grid on; grid minor

% grid size
Dx = 10;
Dz = 10;

% create grids
XX = x_min:Dx:x_max;
ZZ = z_min:Dz:z_max;
LX = length(XX);
LZ = length(ZZ);

MC = zeros(LZ,LX);

set(gca,'YDir','reverse')
%% add first boundary
[x_Bn,z_Bn] = ginput;

flag_fill = 'upper'; % either 'upper' or 'lower'
interp_type = 'pchip'; % linear, spline, nearest, pchip
flag_exlat = 'extd' ; % extd (extend) or ctd (contained)

if strcmp(flag_exlat,'extd')
    XBn = XX; ZBn = interp1(x_Bn,z_Bn,XBn,interp_type);
elseif strcmp(flag_exlat,'ctd')
    XBn = min(x_Bn):Dx:max(x_Bn);
    ZBn = interp1(x_Bn,z_Bn,XBn,interp_type);
end
a_val = find(~isnan(ZBn));
z_Bn = ZBn(a_val);
x_Bn = XBn(a_val);
ZBn = interp1(x_Bn,z_Bn,XBn,interp_type); 

ZBni = round(ZBn/Dz)+1; %convert value to index
ZBni(ZBni<1) = 1;
ZBni(ZBni>LZ) = LZ;

XBni = round(XBn/Dx)+1; %convert value to index
XBni(XBni<x_min) = x_min;
XBni(XBni>x_max) = x_max;

Lke = 0;
% add a horizontal boundary
if strcmp(flag_fill,'upper')
    for ix = 1:length(ZBn)
        MC(1:ZBni(ix),XBni(ix)) = Lke+1;
    end
elseif strcmp(flag_fill,'lower')
    for ix = 1:length(ZBn)
        MC(ZBni(ix):end,XBni(ix)) = Lke+1;
    end
end

imagesc(XX,ZZ,MC); hold on
plot(XX,ZBn,'LineWidth',2)

xlim([x_min x_max])
ylim([z_min z_max])
set(gca,'YDir','reverse')
grid on; grid minor
colormap 'bone'
Lke = Lke+1;

%% add another boundaries
Lke = Lke+1; 
[x_Bn,z_Bn] = ginput;

flag_fill = 'upper'; % 'upper' or 'lower'
interp_type = 'pchip'; % pchip, spline, linear, nearest
flag_exlat = 'extd' ; % extd (extend) or ctd (contained)

if strcmp(flag_exlat,'extd')
    XBn = XX; ZBn = interp1(x_Bn,z_Bn,XBn,interp_type);
elseif strcmp(flag_exlat,'ctd')
    XBn = min(x_Bn):Dx:max(x_Bn);
    ZBn = interp1(x_Bn,z_Bn,XBn,interp_type);
end
a_val = find(~isnan(ZBn));
z_Bn = ZBn(a_val);
x_Bn = XBn(a_val);
ZBn = interp1(x_Bn,z_Bn,XBn,interp_type); 

ZBni = round(ZBn/Dz)+1; %convert value to index
ZBni(ZBni<1) = 1;
ZBni(ZBni>LZ) = LZ;

XBni = round(XBn/Dx)+1; %convert value to index
XBni(XBni<x_min) = x_min;
XBni(XBni>x_max) = x_max;

% add a horizontal boundary
if strcmp(flag_fill,'upper')
    for ix = 1:length(ZBn)
        MC(1:ZBni(ix),XBni(ix)) = Lke;
    end
elseif strcmp(flag_fill,'lower')
    for ix = 1:length(ZBn)
        MC(ZBni(ix):end,XBni(ix)) = Lke+1;
    end
end

imagesc(XX,ZZ,MC); hold on
plot(XBn,ZBn,'LineWidth',2)

xlim([x_min x_max])
ylim([z_min z_max])
set(gca,'YDir','reverse')
grid on; grid minor

%% add a polygon
Lke = Lke+1;
[x_Bn,z_Bn] = ginput;
interp_type = 'linear';
XBn = [];ZBn = [];
for ip = 2:length(x_Bn)
    XBnp = min(x_Bn(ip-1:ip)):Dx:max(x_Bn(ip-1:ip));
    ZBnp = interp1(x_Bn(ip-1:ip),z_Bn(ip-1:ip),XBnp,interp_type);
    XBnp = round(XBnp/Dx);
    ZBnp = round(ZBnp/Dz);
    unikz = find(diff(ZBnp)~=0);
    XBn = [XBn XBnp(unikz)]; 
    ZBn = [ZBn ZBnp(unikz)]; 
    figure(999);
    plot(XBnp,ZBnp,'LineWidth',2); hold on; grid on; grid minor;
    xlim([x_min x_max])
    ylim([z_min z_max])
    set(gca,'YDir','reverse')
end

    XBnp = min(x_Bn(end),x_Bn(1)):Dx:max(x_Bn(end),x_Bn(1));
    ZBnp = interp1([x_Bn(end) x_Bn(1)],[z_Bn(end) z_Bn(1)],XBnp,interp_type);
    XBnp = round(XBnp/Dx);
    ZBnp = round(ZBnp/Dz);
    unikz = find(diff(ZBnp)~=0);
    XBn = [XBn XBnp(unikz)]; 
    ZBn = [ZBn ZBnp(unikz)]; 

    figure(999);
    plot(XBnp,ZBnp,'LineWidth',2); hold on; grid on; grid minor;
    xlim([1 LX])
    ylim([1 LZ])
    set(gca,'YDir','reverse')
figure(999);
plot(XBn,ZBn,'.k'); hold on; grid on; grid minor;
xlim([1 LX])
ylim([1 LZ])
set(gca,'YDir','reverse')

% fill 1
for iz = 1:LZ
    z_ip = find(ZBn==iz);
    if any(z_ip)
        for ix = 1:LX
            if any(XBn(z_ip)>ix) && any(XBn(z_ip)<ix)
                if mod(length(find(XBn(z_ip)>ix)),2)==1 ||...
                        mod(length(find(XBn(z_ip)<ix)),2)==1
                    MC(iz,ix) = Lke+1;
%                     MC(iz,ix) = 21; % use to replace
                end             
            end
        end
    end
end

for iz = 2:LZ-1
    for ix = 2:LX-1
        if MC(iz,ix-1) == MC(iz,ix+1) && MC(iz,ix) ~= MC(iz,ix+1)
            MC(iz,ix-1) = MC(iz,ix);
            MC(iz,ix+1) = MC(iz,ix);
        end
    end
end

figure(7485)
imagesc(XX,ZZ,MC); hold on
set(gca,'YDir','reverse')
grid on; grid minor

%% add a smooth polygon
Lke = Lke+1;
[x_Bn,z_Bn] = ginput;
interp_type = 'spline';
XBn = [];ZBn = [];
for ip = 2:length(x_Bn)-1
    XBnp = min(x_Bn(ip-1:ip)):Dx:max(x_Bn(ip-1:ip));
    ZBnp = interp1(x_Bn(ip-1:ip+1),z_Bn(ip-1:ip+1),XBnp,interp_type);
    XBnp = round(XBnp/Dx);
    ZBnp = round(ZBnp/Dz);
    unikz = find(diff(ZBnp)~=0);
    XBn = [XBn XBnp(unikz)]; 
    ZBn = [ZBn ZBnp(unikz)]; 
    figure(999);
    plot(XBnp,ZBnp,'LineWidth',2); hold on; grid on; grid minor;
    xlim([x_min x_max])
    ylim([z_min z_max])
    set(gca,'YDir','reverse')
end

    XBnp = min(x_Bn(end-1),x_Bn(1)):Dx:max(x_Bn(end-1),x_Bn(1));
    ZBnp = interp1([x_Bn(end-1) x_Bn(end) x_Bn(1)],[z_Bn(end-1) z_Bn(end) z_Bn(1)],XBnp,interp_type);
    XBnp = round(XBnp/Dx);
    ZBnp = round(ZBnp/Dz);
    unikz = find(diff(ZBnp)~=0);
    XBn = [XBn XBnp(unikz)]; 
    ZBn = [ZBn ZBnp(unikz)]; 

    figure(999);
    plot(XBnp,ZBnp,'LineWidth',2); hold on; grid on; grid minor;
    xlim([1 LX])
    ylim([1 LZ])
    set(gca,'YDir','reverse')
figure(999);
plot(XBn,ZBn,'.k'); hold on; grid on; grid minor;
xlim([1 LX])
ylim([1 LZ])
set(gca,'YDir','reverse')

% fill 1
for iz = 1:LZ
    z_ip = find(ZBn==iz);
    if any(z_ip)
        for ix = 1:LX
            if any(XBn(z_ip)>ix) && any(XBn(z_ip)<ix)
                if mod(length(find(XBn(z_ip)>ix)),2)==1 ||...
                        mod(length(find(XBn(z_ip)<ix)),2)==1
                    MC(iz,ix) = Lke+1;
%                     MC(iz,ix) = 12;
                end             
            end
        end
    end
end

for iz = 2:LZ-1
    for ix = 2:LX-1
        if MC(iz,ix-1) == MC(iz,ix+1) && MC(iz,ix) ~= MC(iz,ix+1)
            MC(iz,ix-1) = MC(iz,ix);
            MC(iz,ix+1) = MC(iz,ix);
        end
    end
end

figure(7485)
imagesc(XX,ZZ,MC); hold on
set(gca,'YDir','reverse')
grid on; grid minor

%% Create medium
Param_1 = zeros(size(MC));
Param_2 = zeros(size(MC));

figure(629)
subplot(2,1,1)
imagesc(ZZ,XX,Param_1); grid on; grid minor; colorbar
title 'Param1'
subplot(2,1,2)
imagesc(ZZ,XX,Param_2); grid on; grid minor; colorbar
title 'Param2'
%% Add medium parameter
Lke = 1;
Param_1(MC== Lke ) = 2000;
Param_2(MC== Lke ) = 3200;

figure(629)
subplot(2,1,1)
imagesc(XX,ZZ,Param_1); grid on; grid minor; colorbar
title 'Param1'
subplot(2,1,2)
imagesc(XX,ZZ,Param_2); grid on; grid minor; colorbar
title 'Param2'

