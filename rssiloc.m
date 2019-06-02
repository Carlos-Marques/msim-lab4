%
% Localize a source with unknown power from RSSI measurements
%

% Simulation setup
N = 10;					% No. of anchors
n = 2;					% Embedding dimension
sidelength = 100;
a = sidelength*rand(n,N);		% Anchor positions
x = sidelength*rand(n,1);		% Source position

D = squareform(pdist([x zeros(size(x)) a]'));
d = D(1,3:end);				% Source-anchor distances
an = D(2,3:end);			% Anchor norms

% Generate observations
P0 = 100;				% Source power
P = P0./(d.^2);				% Noiseless RSSI

stdev = 1e-1;				% Log-noise standard deviation
%stdev = 0;
P = P.*exp(stdev*randn(size(P)));	% Introduce noise
QP = 1e-2;
P = QP*round(P/QP);			% Quantize power measurements

% Localize source by least-squares
A = [-2*repmat(P,[n 1]).*a; -ones(size(P)); P]';
b = (-P.*(an.^2))';

z = A\b;
xe = z(1:n);
norm(x-xe)

if n==2,
  plot(a'*[1; 1i],'o'); hold all
  plot(x'*[1; 1i],'x'); plot(xe'*[1; 1i],'s'); hold off
  axis(sidelength*[0 1 0 1]); axis('square')
end

% RLS formulation (one-shot)
RlsPar = struct('lam',1);
[e,w,RlsPar] = qrrls(A,b,RlsPar);
norm(z-w)

% RLS formulation (incremental)
RlsPar = struct('lam',1);
for i = 1:size(A,1)
  [e,w,RlsPar] = qrrls(A(i,:),b(i),RlsPar);
end
norm(z-w)
