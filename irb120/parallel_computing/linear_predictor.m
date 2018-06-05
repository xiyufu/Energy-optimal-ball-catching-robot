function [x_inter, si] = linear_predictor(a_opt,b_opt,dt,si, s, c, fs, optis, optis_param, path_param, path)

    % global f_m f_c f_g f_fc Sign solver;
    f_m = fs{1};
    f_c = fs{2};
    f_g = fs{3};
    f_fc = fs{4};
    Sign = fs{5};
    solver = fs{6};
    % global optis
    % global q_meas dqdt_meas b_optimized si_old si_orig grid_size si_new
    q_meas = optis_param{1};
    dqdt_meas = optis_param{2};
    b_optimized = optis_param{3};
    si_old = optis_param{4};
    si_orig = optis_param{5};
    grid_size = optis_param{6};
    si_new = optis_param{7};
    % global c0op c1op c2op c3op c4op
    c0op = path_param{1};
    c1op = path_param{2};
    c2op = path_param{3};
    c3op = path_param{4};
    c4op = path_param{5};
    % global qs dqs ddqs;
    qs = path{1};
    dqs = path{2};
    ddqs = path{3};
    
    fs = 250;
    Ts = 1/fs;
    number_of_execution = floor(dt(1)/Ts);
    if dt(1) < Ts
        number_of_execution = 1;
    end
    x_inter = zeros(12,number_of_execution+1); % states during internal loop
    x_inter(:,1) = zeros(12,1);%x(:,k); % set initial values
    c0 = c(1);
    c1 = c(2);
    c2 = c(3);
    c3 = c(4);
    c4 = c(5);
    b_two = [b_opt(1),b_opt(2)];
    optis.set_value(b_optimized, b_two ); % Unchanged inside the loop ...
    optis.set_value(c0op, c0); % so we pass these values to optis now
    optis.set_value(c1op, c1);
    optis.set_value(c2op, c2);
    optis.set_value(c3op, c3);
    optis.set_value(c4op, c4);
    ds_k = s(2:end) - s(1:end-1);
    optis.set_value(grid_size, ds_k(1));
    optis.set_value(si_orig, s(1));
    nu = zeros(6,1); % state of integrator in PID control
    slope_b = (b_opt(2) - b_opt(1))/ds_k(1);
    for kkk = 1:number_of_execution
        b_inter = slope_b*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % q, dqdt, ddqdtt
        q_inter = x_inter(1:6,kkk);
        dq_inter = x_inter(7:12,kkk)/sqrt(b_inter);
        ddq_inter = ddqs(si);        
        % reference
        ds_inter = Ts*Ts*slope_b/4 + Ts*sqrt(b_inter);
        b_inter_next = slope_b*ds_inter + b_inter;
        q_ref_inter = qs(si+ds_inter);
        dq_ref_inter = dqs(si+ds_inter)*sqrt(b_inter_next);
        % m c g
        temp_m = full( f_m(q_inter) );
        temp_c = full( f_c(q_inter, dq_inter) );
        temp_g = full( f_g(q_inter, [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        m_inter = temp_m*dq_inter;
        c_inter = temp_m*ddq_inter + temp_c*dq_inter;
        g_inter = temp_g + f_fc*Sign(dq_inter);
        % get tau
        u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
%         saturation_index = u > tau_max;
%         u(saturation_index) = tau_max(saturation_index);
%         saturation_index = u < -tau_max;
%         u(saturation_index) = -tau_max(saturation_index);
        try % simulation
            out = solver('x0',[x_inter(:, kkk);nu],'p',[u;q_ref_inter;dq_ref_inter]);          
        catch
            q_bad = x_inter(1:6,kkk);
            warning('IDAS failed');
            flag_integratorF = 1;
            break;
        end
        % % % Recording control histroy
%         u_act(:, i_act) = u./tau_max;
%         i_act = i_act + 1;
        % % % end recording 
        temp_state = full(out.xf); % update states
        x_inter(:, kkk+1) = temp_state(1:12);
        nu = temp_state(13:end);
        % update si
        optis.set_value(si_old, si);
        optis.set_value(q_meas, x_inter(1:6, kkk+1));
%         optis.set_value(dq_meas, x_inter(7:12, i+1)/sqrt(b_inter_next));
        optis.set_value(dqdt_meas, x_inter(7:12, kkk+1));
        optis.set_value(b_optimized, [b_opt(1),b_opt(2)]);
        optis.set_value(c0op, c0);
        optis.set_value(c1op, c1);
        optis.set_value(c2op, c2);
        optis.set_value(c3op, c3);
        optis.set_value(c4op, c4);
        optis.set_initial(si_new, si);
        sol_s = optis.solve();
        si = sol_s.value(si_new);
    
        warning('LP finished %d', kkk);

    end % end of execution

end

