#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

#include <visualization.h>
#include <init_state.h>
#include <find_min_vertices.h>
#include <find_max_vertices.h>
#include <fixed_point_constraints.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <T_linear_tetrahedron.h>
#include <V_linear_tetrahedron.h>
#include <V_spring_particle_particle.h>
#include <dV_linear_tetrahedron_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <mass_matrix_mesh.h>
#include <assemble_forces.h>
#include <assemble_stiffness.h>
#include <linearly_implicit_euler.h>
#include <implicit_euler.h>
#include <build_skinning_matrix.h>
#include <dirichlet_boundary_condition_stiky.h>
#include <bounding_box.h>
#include <init_elasticity_stiffness.h>
#include <penalty_force_collision.h>
#include <boundary.h>
#include <average_edge_length.h>
#include <testTools.h>
#include <transform.h>

//Variable for geometry
Eigen::MatrixXd V; //vertices of simulation mesh 
Eigen::MatrixXi T; //tetrahedrons of simulation mesh
Eigen::MatrixXi F; //faces of simulation mesh

// for multi body simulation
std::vector<SoftBody> softBodies;
std::vector<int> pointers;
std::vector<Boundary> boundaryList; // boundary of the entire scene such as wall and floor

// for hash configuration
int numX = 1;
int numY = 1;
int numZ = 1;
double gridSpacing = 1;

//variables for skinning
Eigen::MatrixXd V_skin;
Eigen::MatrixXi F_skin;
Eigen::SparseMatrixd N; 

//material parameters
double density = 0.1;
double YM = 6e5; //young's modulus
double mu = 0.4; //poissons ratio
double D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
double C = 0.5*YM/(2.0*(1.0+mu));

//BC (dirichlet boundary conditions)
// sticky
std::vector<size_t> fixed_point_indices;
std::vector<std::vector<size_t>> fixedPointIndexList; // list for multi bodies
Eigen::SparseMatrixd P; // select matrix
Eigen::VectorXd x0; // positions of fixed points

//mass matrix
Eigen::SparseMatrixd M;
//tetrahedron volumes
Eigen::VectorXd v0;
// damping configuration
double damping = 0.001;
Eigen::SparseMatrixd B;
// preserve the value of M^-1*f_ext for computing IP
Eigen::SparseLU<Eigen::SparseMatrixd, Eigen::COLAMDOrdering<int>> LUSolveM;
Eigen::VectorXd Minv_fext;

//temporal variables
Eigen::VectorXd tmp_q;
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::SparseMatrixd tmp_stiffness; // it can be initialized after the geometry has been loaded
Eigen::MatrixXi elasticityHessianMap;
std::vector<int> matrixMap; // for fast regulating the matrix size cause of constraints from Dirichlet boundary conditions
Eigen::SparseMatrixd tmpMat; // for fast regulating the matrix size

// for forces
std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points;
Eigen::VectorXd externalForce;
Eigen::VectorXd dampingForce;
Eigen::VectorXd boundaryForce;

// comtrol paremeters
bool skinning_on = true;
bool fully_implicit = true;
bool bunny = true; 

// for concurrent computing
Eigen::MatrixXi Tk; // preserve ptr of each vertex in a tetrahedron for force clusters
std::vector<std::vector<Eigen::Vector3d>> forceClusters; // for preserving force
std::vector<std::map<int, Eigen::Matrix3d>> hessianClusters; // for preserving Hessian

//selection spring
double k_selected = 1e5;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {

    double V_ele, T_ele, KE, PE;

    spring_points.clear();

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);

    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {
        spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
    }

    // std::cout << spring_points.size() << std::endl;
    if (spring_points.size() > 0) {
        for (size_t i = 0; i < spring_points.size(); i++) {
            std::cout << "the id is " << spring_points[i].second << std::endl;
        }
    }

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> x)->double {
        double E = 0;
        Eigen::VectorXd fullq = P.transpose()*(q+dt*x)+x0;
        Eigen::VectorXd newq = q + dt * x;

        std::chrono::high_resolution_clock::time_point start, end;
        std::chrono::duration<double, std::milli> localTime = end - start;

        // serial accumulate potential energy
        /*
        double testE = 0;
        start = std::chrono::high_resolution_clock::now();
        for(unsigned int ei=0; ei<T.rows(); ++ei) {
            V_linear_tetrahedron(V_ele, fullq, V, T.row(ei), v0(ei), C, D);
            testE += V_ele;
        }
        end = std::chrono::high_resolution_clock::now();
        localTime = end - start;
        std::cout << "time of computing energy of a tet = " << localTime.count() / T.rows() << std::endl;
        std::cout << "time of computing energy of the system = " << localTime.count() << std::endl;
         */

        // concurrent
        struct SumV {
            double sumV;
            const Eigen::VectorXd& x;

            void operator()(const tbb::blocked_range<size_t>& r) {
                double sum = sumV;
                for (size_t i = r.begin(); i < r.end(); i++) {
                    double e = 0.;
                    V_linear_tetrahedron(e, x , V, T.row(i), v0(i), C, D);
                    sum += e;
                }
                sumV = sum;
            }

            SumV(SumV& S, tbb::split)
                    : sumV(0.), x(S.x) {}

            SumV(const Eigen::VectorXd& otherX)
                    : x(otherX), sumV(0.) {}

            void join (const SumV& other) {sumV += other.sumV;}
        };
        SumV sv(fullq);
        start = std::chrono::high_resolution_clock::now();
        tbb::parallel_reduce(tbb::blocked_range<size_t>(0, T.rows()), sv);
        end = std::chrono::high_resolution_clock::now();
        localTime = end - start;
        // std::cout << "time of parallel_reduce = " << localTime.count() << std::endl;
        E = sv.sumV;

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {   
            V_spring_particle_particle(V_ele, spring_points[pickedi].first, fullq.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            E += V_ele;
        }

        // without external force
        // E += 0.5*(x-qdot).transpose()*M*(x-qdot);

        // with external force
        E *= dt * dt;
        Eigen::VectorXd qhat = q + dt * qdot + dt * dt * P * Minv_fext;
        E += 0.5 * (newq - qhat).transpose() * M * (newq - qhat);

        return E;
    };

    auto force = [&](
            Eigen::VectorXd &f,
            Eigen::Ref<const Eigen::VectorXd> q2,
            Eigen::Ref<const Eigen::VectorXd> qdot2) {
        Eigen::VectorXd testForce(f);
        // elastic force
        std::chrono::high_resolution_clock::time_point start, end;
        start = std::chrono::high_resolution_clock::now();
        assembleForcesFast(f, P.transpose()*q2+x0, V, T, v0, forceClusters, Tk, C, D);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> frameTime = end - start;
        // std::cout << "new version time of Vector = " << frameTime.count() << std::endl;
        // start = std::chrono::high_resolution_clock::now();
        // assemble_forces(testForce, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D);
        // end = std::chrono::high_resolution_clock::now();
        // frameTime = end - start;
        // std::cout << "old version time = " << frameTime.count() << std::endl;

        // test force
        /*
        if (f.isApprox(testForce)) {
            std::cout << "the concurrently assembling forces is correct " << std::endl;
        } else {
            std::cout << "the concurrently assembling forces has error " << std::endl;
        }
         */

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
            dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose()*q2+x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
        }

        // add external forces
        f += externalForce;

        f = P*f;
    };

    //assemble stiffness matrix,
    auto stiffness = [&](
            Eigen::SparseMatrixd &K,
            Eigen::Ref<const Eigen::VectorXd> q2,
            Eigen::Ref<const Eigen::VectorXd> qdot2) {
        Eigen::SparseMatrixd testK(K);
        std::chrono::high_resolution_clock::time_point start, end;
        std::chrono::duration<double, std::milli> frameTime = end - start;
        start = std::chrono::high_resolution_clock::now();
        assembleHessianFast(K, P.transpose()*q2+x0, V, T, v0, hessianClusters, C, D);
        end = std::chrono::high_resolution_clock::now();
        frameTime = end - start;
        // std::cout << "concurrent version time of Matrix = " << frameTime.count() << std::endl;
        // start = std::chrono::high_resolution_clock::now();
        // assemble_stiffness(testK, elasticityHessianMap, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D);
        // end = std::chrono::high_resolution_clock::now();
        // frameTime = end - start;
        // std::cout << "serial version time = " << frameTime.count() << std::endl;

        // test assembling stiffness
        /*
        if (testK.isApprox(K)) {
            std::cout << "the concurrently assembling stiffness is correct " << std::endl;
        } else {
            std::cout << "the concurrently assembling stiffness has error " << std::endl;
        }
         */

        K *= dt * dt;

        testK = K;

        start = std::chrono::high_resolution_clock::now();
        // constraintSize(testK, tmpMat, matrixMap);
        /*
        if (K.isApprox(testK)) {
            std::cout << "the concurrent regulating size is correct " << std::endl;
        } else {
            std::cout << "the concurrent regulating size has error " << std::endl;
        }
         */
        K = P*K*P.transpose();
        end = std::chrono::high_resolution_clock::now();
        frameTime = end - start;
        // std::cout << "regulate size time = " << frameTime.count() << std::endl;
    };

    std::cout << "start solver at " << t * 100 << std::endl;
    if(fully_implicit)
        implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_q, tmp_qdot, tmp_force, tmp_stiffness);
    else
        linearly_implicit_euler(q, qdot, dt, M, force, stiffness, tmp_force, tmp_stiffness);
    /*
    KE = 0;
    PE = 0;

    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D);
        PE += V_ele;
    }
    
    Visualize::add_energy(t, KE, PE);
    */

    // test::visualVector(q);
}

void gravitySim(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {
    double V_ele, T_ele, KE, PE;

    // std::cout << "q : " << std::endl;
    // test::visualVector(q);
    // test::visualVectorEle(q, 1);
    // std::cout << "qdot : " << std::endl;
    // test::visualVector(qdot);

    // contact
    CollidersPtr colliders = std::make_shared<PenaltyForceColliders>(softBodies, pointers, boundaryList);
    SpatialHash table(numX, numY, numZ, gridSpacing);

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> x)->double {
        double E = 0;
        double KE = 0;
        Eigen::VectorXd fullq = P.transpose()*(q+dt*x)+x0; // change the parameter from qdot_1 to x
        Eigen::VectorXd newq = q + dt * x;
        // Eigen::VectorXd newq = P.transpose() * x + x0;

        // serial reduce for summing V_ele
        /*
        for(unsigned int ei=0; ei<T.rows(); ++ei) {

            V_linear_tetrahedron(V_ele, newq , V, T.row(ei), v0(ei), C, D);
            E += V_ele;

            // std::cout << "V" << ei << "=" << V_ele << std::endl;
        }
         */

        // concurrent
        struct SumV {
            double sumV;
            const Eigen::VectorXd& x;

            void operator()(const tbb::blocked_range<size_t>& r) {
                double sum = sumV;
                for (size_t i = r.begin(); i < r.end(); i++) {
                    double e = 0.;
                    V_linear_tetrahedron(e, x , V, T.row(i), v0(i), C, D);
                    sum += e;
                }
                sumV = sum;
            }

            SumV(SumV& S, tbb::split)
                    : sumV(0.), x(S.x) {}

            SumV(const Eigen::VectorXd& otherX)
                    : x(otherX), sumV(0.) {}

            void join (const SumV& other) {sumV += other.sumV;}
        };
        SumV sv(fullq);
        // tbb::parallel_reduce(tbb::blocked_range<size_t>(0, T.rows()), sv);
        // E = sv.sumV;

        // serial for multi bodies' elastic energy
        for (int id = 0; id < colliders->objects.size(); id++) {
            const SoftBody& obj = colliders->objects[id];
            Eigen::Ref<const Eigen::MatrixXi> T = obj.T;
            Eigen::Ref<const Eigen::MatrixXd> V = obj.V;
            int initialIndex = colliders->pointers[id];
            Eigen::Ref<const Eigen::VectorXd> localQ = fullq.segment(initialIndex, 3*V.rows());
            for (int ei = 0; ei < T.rows(); ei++) {
                V_linear_tetrahedron(V_ele, localQ, V, T.row(ei), obj.volumes[ei], C, D);
                E += V_ele;
            }
        }
        // std::cout << "ElaE = " << E << std::endl;

        // drag force which decreases the system energy
        Eigen::VectorXd drag = -damping * x;
        E -= newq.transpose() * drag;

        // contact energy
        E += colliders->collisionEnergy(table, newq, P.transpose()*x);
        // std::cout << "collE = " << dt * dt * colliders->collisionEnergy(table, newq, P.transpose()*x) << std::endl;

        E *= dt * dt;

        // for parts of acceleration of external forces from collision between objects and boundaries
        /*
        Eigen::VectorXd accOfContact;
        accOfContact.resize(q.size());
        accOfContact.setZero();
        accOfContact = LUSolveM.solve(boundaryForce);
        Minv_fext += accOfContact;
        */

        // velocity based method
        // KE = 0.5*(x-qdot).transpose()*M*(x-qdot);
        // position based method
        Eigen::VectorXd qhat = q + dt * qdot + dt * dt * P * Minv_fext; // with gravity and other external forces
        // Eigen::VectorXd qhat = q + dt * qdot; // without gravity
        KE = 0.5 * (newq - qhat).transpose() * M * (newq - qhat);

        // std::cout << "KE = " << KE << std::endl;

        E += KE;

        // std::cout << "full E = " << E << std::endl;
        return E;
    };

    auto force = [&](
            Eigen::VectorXd &f,
            Eigen::Ref<const Eigen::VectorXd> q2,
            Eigen::Ref<const Eigen::VectorXd> qdot2) {
        // assemble elastic forces
        // assemble_forces(f, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D);
        assembleMultiForces(f, P.transpose()*q2+x0, colliders, C, D);
        // std::cout << "elastic force: " << std::endl;
        // test::visualVector(f);

        // add external forces
        f += externalForce;
        colliders->boundaryForce(boundaryForce, table, P.transpose()*q2+x0, P.transpose()*qdot2); // contact with boundaries
        f += boundaryForce;
        // test::visualVector(boundaryForce);

        // contact forces
        colliders->collisionForce(f, table, P.transpose()*q2+x0, P.transpose()*qdot2); // contact with other bodies
        // std::cout << "collision force" << std::endl;
        // test::visualVector(f);

        // drag forces
        dampingForce = -damping * P.transpose() * qdot2;
        f += dampingForce;

        // std::cout << "force: " << std::endl;
        // test::visualVector(f);
        f = P*f;
    };

    //assemble stiffness matrix,
    auto stiffness = [&](
            Eigen::SparseMatrixd &K,
            Eigen::Ref<const Eigen::VectorXd> q2,
            Eigen::Ref<const Eigen::VectorXd> qdot2) {
        // assemble_stiffness(K, elasticityHessianMap, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D);
        assembleMultiStiffness(K, P.transpose()*q2+x0, colliders, C, D);
        K *= dt * dt;

        // test::visualSparseMatrix(K);

        // contact stiffness
        Eigen::SparseMatrixd Hq(K.rows(), K.cols());
        Hq.setZero();
        colliders->collisionStiffness(Hq, P.transpose()*q2+x0);
        Hq *= dt * dt;
        K += Hq;

        // contact damping
        Eigen::SparseMatrixd Hv(K.rows(), K.cols());
        Hv.setZero();
        colliders->collisionDamping(Hv, P.transpose()*qdot2);
        Hv *= dt;
        K += Hv;

        // drag forces
        K -= dt * B;

        K = P*K*P.transpose();

        if(test::isMatSymmetric(K)) {
            // std::cout << "symmetric" << std::endl;
        } else {
            std::cout << "complete stiffness is not symmetric" << std::endl;
            // test::visualSparseMatrix(K);
            // exit(0);
        }
    };

    // for collision construction
    auto collision = [&](Eigen::Ref<const Eigen::VectorXd> v) {
        Eigen::VectorXd x = q + dt * v;
        colliders->initialize();
        colliders->computeConstraintSet(table, P.transpose()*x+x0);
    };

    std::cout << "start solver " << t*100 << "----------------" << std::endl;
    // test::visualVector(q);
    // test::visualVector(qdot);
    if(fully_implicit) {
        // implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_q, tmp_qdot, tmp_force, tmp_stiffness);
        implicitEuler(q, qdot, dt, M, energy, force, stiffness, collision, tmp_qdot, tmp_force, tmp_stiffness);
    }

    /*
    KE = 0;
    PE = 0;
    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D);
        PE += V_ele;
    }

    Visualize::add_energy(t, KE, PE);
     */
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {
    // update vertex positions using simulation for single body
    // Visualize::update_vertex_positions(0, P.transpose() * q + x0);

    // update vertex positions using simulation for multi bodies
    Eigen::VectorXd newq = P.transpose() * q + x0;
    for (int id = 0; id < softBodies.size(); id++) {
        int initialIndex = pointers[id];
        Eigen::Ref<const Eigen::VectorXd> p = newq.segment(initialIndex, 3*softBodies[id].V.rows());
        Visualize::update_vertex_positions(id, p);
    }
}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    if(key =='N') {
        std::cout<<"toggle integrators \n";
        fully_implicit = !fully_implicit;
    } else if(key == 'S') {
        
        skinning_on = !skinning_on;
        Visualize::toggle_skinning(skinning_on);
    }

    return false;
}

inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 
    igl::readMESH("../data/coarser_bunny.mesh",V,T, F);
    igl::readOBJ("../data/bunny_skin.obj", V_skin, F_skin);

    if(argc > 1) {
        if(strcmp(argv[1], "arma") == 0) {
            read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");
            igl::readOBJ("../data/armadillo.obj", V_skin, F_skin);
        
            bunny = false;
            fully_implicit = true;
        }
    }
    
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();
    
    build_skinning_matrix(N, V, T, V_skin);

    //setup simulation 
    init_state(q,qdot,V);

    //add geometry to scene
    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::toggle_skinning(false);
    
    //bunny
    if(bunny)
        Visualize::set_picking_tolerance(1.);
    else
        Visualize::set_picking_tolerance(0.01);

    //volumes of all elements
    igl::volume(V,T, v0);

    //Mass Matrix
    mass_matrix_mesh(M, qdot, T, density, v0);

    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }
    
    //setup constraint matrix
    if(bunny)
        find_min_vertices(fixed_point_indices, V, 3);
    else
        find_min_vertices(fixed_point_indices, V, 0.1);

    //material properties
    //bunny
    if(bunny) {
        YM = 6e6; //young's modulus(k)
        mu = 0.4; //poissons ratio(v)
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu)); // λ in sig coarse
        C = 0.5*YM/(2.0*(1.0+mu)); // μ in sig coarse
        k_selected = 1e8;
    } else {
        //arma
        YM = 6e5; //young's modulus
        mu = 0.4; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e5;
    }

    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixed_point_constraints(P, q.rows(), fixed_point_indices);
    
    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else    
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();

    //igl additional menu setup
    // Add content to the default menu window
    /*
    Visualize::viewer_menu().callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * Visualize::viewer_menu().menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
            "Energy Plot", nullptr,
            ImGuiWindowFlags_NoSavedSettings

        );

        ImVec2 min = ImGui::GetWindowContentRegionMin();
        ImVec2 max = ImGui::GetWindowContentRegionMax();

        max.x = ( max.x - min.x ) / 2;
        max.y -= min.y + ImGui::GetItemsLineHeightWithSpacing() * 3;

        Visualize::plot_energy("T", 1, ImVec2(-15,10), ImVec2(0,2e8), ImGui::GetColorU32(ImGuiCol_PlotLines));
        Visualize::plot_energy("V", 2, ImVec2(-15,10), ImVec2(0,2e7), ImGui::GetColorU32(ImGuiCol_HeaderActive));
        Visualize::plot_energy("T+V", 3, ImVec2(-15,10), ImVec2(0,4e8), ImGui::GetColorU32(ImGuiCol_ColumnActive));

        ImGui::End();
    };

    Visualize::viewer().callback_key_down = key_down_callback;
    */
}

#endif

