% Results file from /casadi_with_cpp/autonaut/solver_gen.cpp
% Generated Apr 14 2023 at 12:10:42
for i = 1:7
    chi_d = eval(['chi_d' num2str(i)]);
    optims = eval(['optims' num2str(i)]);
    plot_prediction(optims, chi_d);
end

function plot_prediction(optims, chi_d)
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
    cmd.beta = rad2deg(atan(cmd.v ./ cmd.u));
    cmd.chi = cmd.psi + cmd.beta;
    
    figure(1); clf; 
    hold on; grid on;
    plot(cmd.psi);
    plot(cmd.chi);
    yline(chi_d);
    legend('\psi', '\chi');
    
    xlabel('time interval'); ylabel('units');
    title('state evolution with time')
end