% dt = 1e-3;
% td_s = rand(101,1);
% td_s = zeros(219,1); td_s(round(end/2)) = 10;
% td_s = sin(2*pi*100*[0:dt:dt*(219-1)]);

L = length(td_s);
t_s = 0:dt:dt*(L-1);

if mod(L,2)==0
    td_s(L+1) = td_s(end);
end

L_f = 2^16;
fd_s = fftshift(fft(td_s,L_f)); % fft and shift f0 to the center

fnyq = 1/2/dt;
f_s = linspace(-fnyq,fnyq,L_f);

if mod(L,2)==0
    td_s = td_s(1:end-1);
end

    figure(744);
    subplot(2,1,1);plot(t_s,td_s);xlabel 'time(s)'; grid on; title 'time domain signal'
    xlim([t_s(1) t_s(end)]); ylim([-1.1*max(abs(td_s)) 1.1*max(abs(td_s))])
    subplot(2,1,2);plot(f_s,abs(fd_s));xlabel 'Hz'; ylabel '|A|'; grid on; title 'Fourier Spectra'
    xlim([f_s(1) f_s(end)]); ylim([0 1.1*max(abs(fd_s))])

