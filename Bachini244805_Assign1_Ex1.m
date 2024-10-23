clc; clearvars; close all;

%% Utility functions

mu = 0.012150;
OM = @(x,y,z, r1, r2)  1/2*(x^2 +y^2)+ (1-mu)/r1 + mu/r2 + 1/2*mu*(1-mu);
Jacobian = @(x, r1, r2) 2*OM(x(1),x(2),x(3), r1, r2) - (x(4)^2+x(5)^2+x(6)^2);
radius = @(L, P) norm(P - L);

%% Ex.1.1

% L4 and L5 from geometrical construction
L4 = [0.5-mu; sqrt(3)/2];
L5 = [0.5-mu; -sqrt(3)/2];

% L1, L2, L3 
dUdx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x+mu-1)./abs(x+mu-1).^3; 
L3 = [fzero(dUdx, -1); 0];
L1 = [fzero(dUdx, 0.5); 0];
L2 = [fzero(dUdx, 1.5); 0];

L_points = [L1 L2 L3 L4 L5];

P1 = [-mu; 0];
P2 = [1-mu; 0];

C_vals = zeros(1, 5);

for ii = 1:5
    L_state = [L_points(:,ii); zeros(4, 1)];
    r1 = radius(L_points(1:2, ii), P1);
    r2 = radius(L_points(1:2, ii), P2);
    C_vals(ii) = Jacobian(L_state, r1, r2);
end

figure
grid on; hold on;
scatter(P1(1), P1(2), 'filled');
scatter(P2(1), P2(2), 'filled');
scatter(L_points(1, :), L_points(2, :), 'filled');
text(L_points(1,:)-0.05, L_points(2,:)-0.075, "L" + (1:5)')
xlim([-1.5 1.5]); ylim([-1 1])

%% Ex 1.2

