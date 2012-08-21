clear all

loaded = load('e_vs_t.dat');
t = loaded(:,1);
E = loaded(:,2);
nw = loaded(:,3);

E_x = 3;
exact = ones(size(E,1),1)*E_x;

%figure(1)
%plot(e_t)

figure(2)
plot(t,E)
hold on
plot(t,exact, 'r')
legend('DMC', 'Exact');
xlabel('\tau')
ylabel('E')

Emin = min(E);
imin = find(E==Emin);
Emax = 2*E_x- Emin;
axis([t(1),t(end), Emin, Emax])
title('2 electron non-interacting case. \Delta\tau = 0.001.')

figure(3)
plot(nw)