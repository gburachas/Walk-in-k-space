function k = archiSpiral(t, slrate, gamma)

A=1;

if nargin == 1,
  slrate = 1;
  gamma = 0.01; %0.01
end;

a = gamma * slrate/A;
b = sqrt(2*a/3);

theta = a* t.^2 ./(2 + (b*t).^(4/3) + 2* (b*t).^(1/3));
r = 1:length(theta);

k = r.* exp(i* theta);

figure(3); clf;
plot(k'); axis equal