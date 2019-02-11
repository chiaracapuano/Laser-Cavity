function [ Hout ] = H(Gamma,cosdt,cavity_length,cavity_loss,var,q,m,n)
% % cavity transfer function approximation in eq.11
c=physconst('LightSpeed')*100;

RTT = 2*cavity_length/c; % round trip time
const = -2/(log(1-cavity_loss)*(2-cavity_loss));
Hout = const*Gamma^2*((1-(1-cavity_loss)*exp(1i*(var-(pi*c)/cavity_length*((q)+2/pi*cosdt*(abs(m)+abs(n)+1)))*RTT))./ ...
    (Gamma^2+(var-(pi*c)/cavity_length*((q)+2/pi*cosdt*(abs(m)+abs(n)+1))).^2));


end

