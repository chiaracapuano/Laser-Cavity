function [ gaussian_pulse, gaussian_pulse_ft  ] = GaussianPulse(A,FWHM,t0,n0,dt,t)

    gamma = 4.0*log(2.0)/FWHM^2;
    norm = A*sqrt(sqrt(gamma/pi));
  
    gaussian_pulse = norm * exp(-0.5*gamma*(t-t0).^2+ 2*pi*1i*(n0)*t);%
    %     plot(gaussian_pulse.^gaussian_pulse)
    gaussian_pulse_ft = fftshift(fft(ifftshift(gaussian_pulse)))*dt;
end

