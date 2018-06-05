function [a_opt,b_opt,dt_opt, s] =NLP_solver(si, T_itc, q_itc, b_init_k, N0, path, fs, opti_vari, opti_param)

    
    % global f_m f_c f_g f_fc Sign;
    f_m = fs{1};
    f_c = fs{2};
    f_g = fs{3};
    f_fc = fs{4};
    Sign = fs{5};
    % global qs dqs ddqs;
    qs = path{1};
    dqs = path{2};
    ddqs = path{3};

%     global opti;
%     global a b tau;
%     global bk_inv;
    opti = opti_vari{1};
    a = opti_vari{2};
    b = opti_vari{3};
    tau = opti_vari{4};
    bk_inv = opti_vari{5};
    % get the collocation points 
    s = linspace(si ,1, N0+1);
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    qh = qs(sh);
    dqh = dqs(sh);
    ddqh = ddqs(sh);
    
    % Evaluate m,c,g at current s 
    m_k = zeros(size(qh));
    c_k = zeros(size(qh));
    g_k = zeros(size(qh));
    for i = 1:N0
        temp_m = full( f_m(qh(:, i)) );
        temp_c = full( f_c(qh(:, i), dqh(:, i)) );
        temp_g = full( f_g(qh(:, i), [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        
        m_k(:, i) = temp_m*dqh(:, i);
        c_k(:, i) = temp_m*ddqh(:, i) + temp_c*dqh(:,i);
        g_k(:, i) = temp_g + f_fc*Sign(dqh(:, i));
    end
    
    tau_holding_k = full( f_g(q_itc, zeros(6,1), zeros(6,1)) );
    % feed data to the solver
    gcp = ones(1,N0);
    
    % global grid_control m c g tau_holding b_init T_remain ds;
    grid_control = opti_param{1};
    m = opti_param{2};
    c = opti_param{3};
    g = opti_param{4};
    tau_holding = opti_param{5};
    b_init = opti_param{6};
    T_remain = opti_param{7};
    ds = opti_param{8};
    
    opti.set_value(grid_control, gcp);
    opti.set_value(m, m_k);
    opti.set_value(c, c_k);
    opti.set_value(g, g_k);
    opti.set_value(tau_holding, tau_holding_k);
    opti.set_value(b_init, b_init_k);
    opti.set_value(T_remain, T_itc);
    opti.set_value(ds, ds_k(1));
    sol = opti.solve();   % actual solve
    
    warning('NLP finished');
    % get the optimal solution
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    tau_opt = sol.value(tau);
    dt_opt = ds_k.*bk_inv_opt;
    q_ref_k= qs(si+ds_k(1));
    dq_ref_k = dqs(si+ds_k(1))*sqrt(b_opt(2));
    
%     t = zeros(1,N0+1);
%     for i = 1:N0
%         t(i+1) = t(i) + dt_opt(i);
%     end
%     plot(t(1:end-1), tau_opt(1,:));
end

