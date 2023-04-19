% Results file from /casadi_with_cpp/autonaut/solver_gen.cpp
% Generated Apr 14 2023 at 12:10:42

nx = 4; nu = 1;
N = 90;
size(optims);

index_states  = nx * (N+1);
index_control = index_states + nu * N;

x_pred = reshape(optims(1:index_states)', nx, N+1)';    
u_pred = reshape(optims(index_states+1:index_control)', nu, N)';           % get controls only from the solution

cmd = struct();
cmd.psi = rad2deg(x_pred(:,1));
cmd.u = x_pred(:,2);
cmd.v = x_pred(:,3);
cmd.r = rad2deg(x_pred(:,4));
cmd.beta = rad2deg(atan(cmd.v./cmd.u));
cmd.chi = cmd.psi + cmd.v;

figure(1); clf; 
hold on; grid on;
plot(cmd.psi);
plot(cmd.chi);
legend('\psi', '\chi');

xlabel('time interval'); ylabel('units');
title('state evolution with time')
