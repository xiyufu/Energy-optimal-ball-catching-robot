// Simple mass example
/*  min  u^2
    s.t. xdd = u/m
         u*xd <= 20
         x(0) == [0, 0]
         x(end) == [10, 0]
         -10 <= u <= 10 
    */
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "acados/utils/print.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_nlp/ocp_nlp_constraints.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/sim/sim_erk_integrator.h"

#include "acados_c/ocp_nlp_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "mass_model/dyn_mass.h"
#include "mass_model/power_limit.h"

int main() {

    int num_states = 2;
    int num_controls = 1;
    int N = 100;
    double Tf = 100.0;
    double R = 1.0;

    int max_num_sqp_iterations = 100;

    double x0[] = {0, 0};
    double x_end[] = {10, 0};
    std::vector<int> nx(N+1, num_states), nu(N+1, num_controls), nbx(N+1, 0), nbu(N+1, num_controls),
                    nb(N+1, num_controls), ng(N+1, 0), nh(N+1, 1), np(N+1, 0),
                    ns(N+1, 0), nv(N+1, num_states+num_controls), ny(N+1, num_states+num_controls);
    //ny is the number of outputs and nv is the number of variables [u;v]
    //in the pendulum example, the output y is v.
    // nh = 1, we have one nonlinear constraints p=u*v < pmax
    nbx.at(0) = num_states; // x(0) = x0
    nb.at(0) = num_states+num_controls; //x(0) = x0; Fmin<u<Fmax
    nb.at(N) = num_states; //x(N) = x_end
    nu.at(N) = 0; // no force at the end
    nv.at(N) = num_states; // nv = nx + nu
    ny.at(N) = num_states; // ny = nv. 

    // make plan
    ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(N);
    plan->nlp_solver = SQP_GN;
    plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
    for (int i = 0; i < N; i++){
        plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        plan->sim_solver_plan[i].sim_solver = ERK;
    }

    ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, N);

    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
    ocp_nlp_dims_initialize(config, nx.data(), nu.data(), ny.data(), nbx.data(), nbu.data(), ng.data(), nh.data(), np.data(), ns.data(), dims);

    external_function_casadi dynamic_model_casadi[N];
    for (int i = 0; i < N; i++){
        dynamic_model_casadi[i].casadi_fun = &dyn_mass;
        dynamic_model_casadi[i].casadi_n_in = &dyn_mass_n_in;
        dynamic_model_casadi[i].casadi_n_out = &dyn_mass_n_out;
        dynamic_model_casadi[i].casadi_sparsity_in = &dyn_mass_sparsity_in;
        dynamic_model_casadi[i].casadi_sparsity_out = &dyn_mass_sparsity_out;
        dynamic_model_casadi[i].casadi_work = &dyn_mass_work;
    }

    // NLP model
    int function_size = 0;
    for (int i = 0; i < N; i++){
        function_size += external_function_casadi_calculate_size(dynamic_model_casadi+i);
    }

    char *c_ptr = (char*) calloc(1, function_size);
    for (int i = 0; i < N; i++){
        external_function_casadi_assign(dynamic_model_casadi+i, c_ptr);
        c_ptr += external_function_casadi_calculate_size(dynamic_model_casadi+i);
    }

    ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

    for (int i = 0; i < N; i++){
        nlp_in->Ts[i] = Tf/N;
    }

    // NLP cost
    ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;

    // cost = 0.5*[u;x]'Cyt*[u;x]
    for (int i = 0; i < N; i++){
        blasfeo_dgese(nv[i], ny[i], 0.0, &cost_ls[i]->Cyt, 0, 0);
        for (int j = 0; j < N; j++){
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, j, nx[i]+j) = 1.0;
        }
        for (int j = 0; j< nx[i]; j++){
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, nu[i]+j, j) = 1.0;
        }

    }

    //W weighting factor of stage cost?
    for (int i = 0; i < N; i++){
        blasfeo_dgese(ny[i], ny[i], 0.0, &cost_ls[i]->W, 0, 0);
        for (int j = 0; j < N; j++){
            BLASFEO_DMATEL(&cost_ls[i]->W, j, j) = 0.0;
            //Q is 0
        }
        for (int j = 0; j < N; j++){
            BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = 1.0;
            //R is 1, cost = u^2
        }
    }

    // WN weigthing factor of terminal cost?
    blasfeo_dgese(ny[N], ny[N], 0.0, &cost_ls[N]->W, 0, 0);
    for (int j = 0; j < nx[N]; j++){
        BLASFEO_DMATEL(&cost_ls[N]->W, j, j) = 10.0;
    }

    // y_ref
    for (int i = 0; i <= N; i++){
        blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);
        //No reference, weighting factor set to 0, perhaps
    }

    // NLP dynamics
    for (int i = 0; i < N; i++){
        ocp_nlp_dynamics_cont_model *dynamics = (ocp_nlp_dynamics_cont_model *) nlp_in->dynamics[i];
        erk_model *model = (erk_model *) dynamics->sim_model;
        model->expl_ode_fun = (external_function_generic *) &dynamic_model_casadi[i];
        // use explict Runge Kutta for shooting method
    }

    nlp_in->freezeSens = false;

    // nonlinear constraints power = u*v <= 20 [W]
    external_function_casadi power_constraint;
    power_constraint.casadi_fun = &plimt;
    power_constraint.casadi_n_in = &plimt_n_in;
    power_constraint.casadi_n_out = &plimt_n_out;
    power_constraint.casadi_sparsity_in = &plimt_sparsity_in;
    power_constraint.casadi_sparsity_out = &plimt_sparsity_out;
    power_constraint.casadi_work = &plimt_work;

    int constraint_size = external_function_casadi_calculate_size(&power_constraint);
    void *ptr = malloc(constraint_size);
    external_function_casadi_assign(&power_constraint, ptr);

    // bounds
    ocp_nlp_constraints_model **states_bounds = (ocp_nlp_constraints_model **) nlp_in->constraints;
    // At initial state: x == x0; Fmin <= u <= Fmax; u*v - 20 <= 0
    double Fmax = 10;
    double Fmin = -Fmax;
    double lb0[] = {Fmin, x0[0], x0[1]};
    double ub0[] = {Fmax, x0[0], x0[1]};
    double lhb = -100000;
    double uhb = 0;
    // x0 <= x(0); Fmin <= u. u and x in what order? u first I guess
    blasfeo_pack_dvec(nb[0], lb0, &states_bounds[0]->d, 0);
    // -100000 <= u*v - pmax
    blasfeo_pack_dvec(nh[0], &lhb, &states_bounds[0]->d, nb[0]+ng[0]);
    // x(0) <= x0; u <= Fmax
    blasfeo_pack_dvec(nb[0], ub0, &states_bounds[0]->d, nb[0]+ng[0]+nh[0]);
    // u*v - pmax <= 0
    blasfeo_pack_dvec(nh[0], &uhb, &states_bounds[0]->d, 2*(nb[0]+ng[0])+nh[0]);
    // push nonlinear constraints to states_bounds
    states_bounds[0]->h = (external_function_generic *) &power_constraint;

    // In the middle, Fmin <= u <= Fmax; -100000 <= u*v - pmax <=0
    for (int i = 1; i<N; i++){
        blasfeo_pack_dvec(nb[i], &Fmin, &states_bounds[i]->d, 0);
        blasfeo_pack_dvec(nh[i], &lhb, &states_bounds[i]->d, nb[i]+ng[i]);
        blasfeo_pack_dvec(nb[i], &Fmax, &states_bounds[i]->d, nb[i]+ng[i]+nh[i]);
        blasfeo_pack_dvec(nh[i], &uhb, &states_bounds[i]->d, 2*(nb[i]+ng[i])+nh[i]);
        states_bounds[i]->h = (external_function_generic *) &power_constraint;
    }

    // At the end, x_end <= x(N) <= x_end; -100000 <= u*v - pmax <= 0
    blasfeo_pack_dvec(nb[N], x_end, &states_bounds[N]->d, 0);
    blasfeo_pack_dvec(nh[N], &lhb, &states_bounds[N]->d, nb[N]+ng[N]);
    blasfeo_pack_dvec(nb[N], x_end, &states_bounds[N]->d, nb[N]+ng[N]+nh[N]);
    blasfeo_pack_dvec(nh[N], &uhb, &states_bounds[N]->d, 2*(nb[N]+ng[N])+nh[N]);
    states_bounds[N]->h = (external_function_generic *) &power_constraint;

    void *nlp_opts = ocp_nlp_opts_create(config, dims);

    ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;
    sqp_opts->maxIter = max_num_sqp_iterations;
    sqp_opts->min_res_g = 1e-9;
    sqp_opts->min_res_b = 1e-9;
    sqp_opts->min_res_d = 1e-9;
    sqp_opts->min_res_m = 1e-9;
    ((ocp_qp_partial_condensing_solver_opts *) sqp_opts->qp_solver_opts)->pcond_opts->N2 = N;

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
    for (int i = 0; i <= N; i++){
        // set initial value to 0
        blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);
    }

    ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    acados_timer timer;
    acados_tic(&timer);

    int solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);

    double elapsed_time = acados_toc(&timer);

    printf("\nsolution\n");
    ocp_nlp_out_print(dims, nlp_out);

    return solver_status;
}



