#include <iostream>
#include <thread>
#include <chrono>

#include <assignment_setup.h>
#include <visualization.h>
#include <oneapi/tbb.h>

#include <testSolver.h>

//simulation state
Eigen::VectorXd q;
Eigen::VectorXd qdot;

//simulation time and time step
double t = 0; //simulation time 
double dt = 0.01; //time step

//simulation loop
bool simulating = true;

bool simulation_callback() {

    while(simulating) {
        std::chrono::high_resolution_clock::time_point start, end;
        start = std::chrono::high_resolution_clock::now();

        // simulate(q, qdot, dt, t); // for a single body
        gravitySim(q, qdot, dt, t); // for multi bodies

        // other module tests
        // testColumnList();
        // testEigenSparse();
        // testSolver();

        t += dt;

        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> frameTime = end - start;
        double frameRate = 1000.0 / frameTime.count();
        // std::cout << "the FPS of an optimization is " << frameRate << std::endl;

        if (t > 0.) {
            // exit(0);
        }
    }

    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    
    draw(q, qdot, t);

    return false;
}

void initialize(Eigen::VectorXd &q, Eigen::VectorXd &qdot) {
    /*
    // load tet
    V.resize(4, 3);
    V << 1.0, 0.0, 0.0,
         0.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    T.resize(1, 4);
    T << 0, 1, 3, 2;
    F.resize(4, 3);
    F << 0, 3, 1,
         0, 2, 3,
         0, 1, 2,
         1, 3, 2;
    */

    // load cube divided by 5 tetrahedrons
    /*
    V.resize(8, 3);
    T.resize(5, 4);
    V << 100.0, 100.0, 0.0,
         100.0, 100.0, 100.0,
         0.0, 100.0, 0.0,
         0.0, 100.0, 100.0,
         100.0, 0.0, 0.0,
         100.0, 0.0, 100.0,
         0.0, 0.0, 100.0,
         0.0, 0.0, 0.0;
    T << 0, 1, 2, 4,
         1, 2, 3, 6,
         1, 4, 5, 6,
         1, 2, 4, 6,
         2, 4, 6, 7;
    */

    /*
    // load cube divided by 6 tetrahedrons
    V.resize(8, 3);
    T.resize(6, 4);
    F.resize(12, 3);
    V << 1.0, 1.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 1.0, 1.0,
         1.0, 1.0, 1.0,
         1.0, 0.0, 0.0,
         0.0, 0.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, 0.0, 1.0;
    T << 0, 3, 1, 5,
         2, 1, 3, 5,
         0, 4, 3, 5,
         3, 5, 4, 7,
         3, 6, 5, 7,
         3, 2, 5, 6;
    F << 2, 3, 1,
         0, 1, 3,
         1, 5, 2,
         6, 2, 5,
         7, 4, 3,
         0, 3, 4,
         2, 6, 3,
         7, 3, 6,
         1, 0, 5,
         4, 5, 0,
         6, 5, 7,
         4, 7, 5;
    */

    /*
    // load a beam constructed by 2 cubes
    V.resize(12, 3);
    T.resize(12, 4);
    V << 0.0, 0.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, 0.0, 1.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 1.0, 1.0,
         1.0, 1.0, 1.0,
         1.0, 1.0, 0.0,
         0.0, 2.0, 0.0,
         0.0, 2.0, 1.0,
         1.0, 2.0, 1.0,
         1.0, 2.0, 0.0;
    T << 0, 4, 6, 7,
         0, 4, 5, 6,
         0, 1, 5, 6,
         0, 1, 2, 6,
         0, 3, 6, 7,
         0, 2, 3, 6,
         4, 8, 10, 11,
         4, 8, 9, 10,
         4, 5, 9, 10,
         4, 5, 6, 10,
         4, 7, 10, 11,
         4, 6, 7, 10;
    */

    // bunny
    // igl::readMESH("../data/coarser_bunny.mesh",V,T, F);
    // for material settings
    {
        density = 5;
        YM = 6e2;
        mu = 0.4;
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e8;
    }

    // armadillo
    // read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");

    // bar
    read_tetgen(V, T, "../data/bar/bar_h.node", "../data/bar/bar_h.ele");

    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();

    // edge set and scale model
    // Eigen::MatrixXi E;
    // getEdgeSetForVolume(E, T);
    // double edge = longestEdge(V, E);
    // double factor = 1.0 / edge;
    // V *= factor;
    BoundingBoxD box(V);
    double factor1 = box.scaleLongestEdge(1);
    V *= factor1;

    // test scale model
    // std::cout << "longest edge is " << edge << std::endl;
    std::cout << "factor is " << factor1 << std::endl;
    // std::cout << box.height() << ", " << box.width() << ", " << box.depth() << std::endl;
    BoundingBoxD newBox(V);
    std::cout << newBox.height() << ", " << newBox.width() << ", " << newBox.depth() << std::endl;

    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(72,209,204)/255.);
    std::cout << "finish add objs to the scene" << std::endl;

    // Visualize::set_picking_tolerance(1.); // for big model
    Visualize::set_picking_tolerance(0.01); // for small model

    init_state(q, qdot, V);

    igl::volume(V, T, v0);
    /*
    for (size_t i = 0; i < v0.size(); i++) {
        // std::cout << "volume" << i << "=" << v0[i] << std::endl;
        if (v0[i] < 0) {
            v0[i] = -v0[i];
        }
    }
     */

    mass_matrix_mesh(M, qdot, T, density, v0);
    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }

    // prepare data for concurrent computing
    // initialize sparse matrix for filling elasticity stiffness
    initElasticityStiffness(tmp_stiffness, elasticityHessianMap, T, q.size());
    // initialize force clusters for concurrently assembling elastic forces
    forceClusters.resize(V.rows());
    hessianClusters.resize(V.rows()*V.rows());
    Tk.resize(T.rows(), 4);
    for (int k = 0; k < T.rows(); k++) {
        for (int i = 0; i < 4; i++) {
            int vi = T(k, i);
            forceClusters[vi].emplace_back(0., 0., 0.);
            Tk(k, i) = forceClusters[vi].size() - 1;

            for (int j = 0; j < 4; j++) {
                int vj = T(k, j);
                int index = vi * V.rows() + vj;
                hessianClusters[index].emplace(k, Eigen::Matrix3d::Zero());
            }
        }
    }

    // set external forces
    // just add gravity to the external forces
    Eigen::Vector3d g(0, -10, 0);
    externalForce.resize(q.size());
    externalForce.setZero();
    for (int i = 0; i < q.size()/3; i++) {
        externalForce.segment(3*i, 3) += g;
    }
    externalForce = M * externalForce;

    // M_inv * f_ext, in fact this is just the total acceleration of external forces
    Eigen::SparseLU<Eigen::SparseMatrixd, Eigen::COLAMDOrdering<int>> LUSolver;
    LUSolver.analyzePattern(M);
    LUSolver.factorize(M);
    Minv_fext.resize(externalForce.size());
    Minv_fext = LUSolver.solve(externalForce);

    std::vector<size_t> bcNodes;
    // readFixedVertices(bcNodes, "../data/bar/fixed_points.txt");
    // findTopBackVertices(bcNodes, V);
    findBackVertices(bcNodes, V, 1e-2);

    // find_min_vertices(fixed_point_indices, V, 3); // for big model
    // find_min_vertices(fixed_point_indices, V, 0.1); // for small model
    // find_max_vertices(fixed_point_indices, V, 0.1);

    // prepare for fast regulating the size of matrix
    // initElasticityStiffness(tmpMat, matrixMap, bcNodes, T, q.size());

    P.resize(q.rows(),q.rows());
    P.setIdentity();
    // fixed_point_constraints(P, q.rows(), fixed_point_indices);
    fixed_point_constraints(P, q.rows(), bcNodes);

    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();
    std::cout << "finish regulate to the right size" << std::endl;

    // Visualize::viewer().callback_key_down = &keydown_callback;
}

void initializeMultiBodySimulation(Eigen::VectorXd &q, Eigen::VectorXd &qdot) {
    // load geometry configurations
    /*
    // load a single tetrahedron
    V.resize(4, 3);
    V << 1.0, 0.0, 0.0,
         0.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    T.resize(1, 4);
    T << 0, 1, 3, 2;
    F.resize(4, 3);
    F << 0, 3, 1,
         0, 2, 3,
         0, 1, 2,
         1, 3, 2;
    */

    /*
    // load cube divided by 6 tetrahedrons
    V.resize(8, 3);
    T.resize(6, 4);
    F.resize(12, 3);
    V << 1.0, 1.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 1.0, 1.0,
         1.0, 1.0, 1.0,
         1.0, 0.0, 0.0,
         0.0, 0.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, 0.0, 1.0;
    T << 0, 3, 1, 5,
         2, 1, 3, 5,
         0, 4, 3, 5,
         3, 5, 4, 7,
         3, 6, 5, 7,
         3, 2, 5, 6;

     F << 2, 3, 1,
         0, 1, 3,
         1, 5, 2,
         6, 2, 5,
         7, 4, 3,
         0, 3, 4,
         2, 6, 3,
         7, 3, 6,
         1, 0, 5,
         4, 5, 0,
         6, 5, 7,
         4, 7, 5;
    */

    // bunny
    // igl::readMESH("../data/coarser_bunny.mesh",V,T, F);

    // armadillo
    read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");

    // bar
    // read_tetgen(V, T, "../data/bar/bar_h.node", "../data/bar/bar_h.ele");

    // get boundary face elements with theoretical outward normal
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();

    // set material parameters
    {
        density = 5;
        YM = 6e2;
        mu = 0.4;
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e8;
    }
    // scale a single object‘s setting to unified size
    BoundingBoxD box(V);
    double factor = box.scaleLongestEdge(1);
    V *= factor;

    SoftBody obj0(V, F, T);
    SoftBody obj1(V, F, T);
    SoftBody obj2(V, F, T);
    obj0.getEdgeSetForVolume();
    obj1.getEdgeSetForVolume();
    obj2.getEdgeSetForVolume();
    Eigen::Vector3d translationUp(0., 1.1, -0.2);
    obj1.translate(translationUp);
    Eigen::Vector3d translationRight(1.1, 0., 0.);
    obj2.translate(translationRight);

    softBodies.push_back(obj0);
    // softBodies.push_back(obj1);
    // softBodies.push_back(obj2);

    pointers.resize(softBodies.size());
    pointers[0] = 0;
    for (size_t i = 1; i < softBodies.size(); i++) {
        pointers[i] = pointers[i-1] + 3 * softBodies[i-1].V.rows();
    }

    V_skin.resize(V.rows(), V.cols());
    F_skin.resize(F.rows(), F.cols());

    for (auto &obj : softBodies) {
        Visualize::add_object_to_scene(obj.V, obj.F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    }
    std::cout << "finish add objs to the scene" << std::endl;

    // for user interaction
    // Visualize::set_picking_tolerance(1.); // for big objects
    // Visualize::set_picking_tolerance(0.01); // for small(unified) objects

    std::vector<Eigen::VectorXd> pos;
    std::vector<Eigen::VectorXd> vel;
    pos.resize(softBodies.size());
    vel.resize(softBodies.size());
    for (size_t i = 0; i < softBodies.size(); i++) {
        softBodies[i].initSelfState(pos[i], vel[i]);
    }
    initMultiState(q, qdot, pos, vel);

    massMatrixMultiBodies(M, qdot, softBodies, pointers, density);
    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }

    // parallel assemble stiffness
    initMultiBodiesStiffness(tmp_stiffness, softBodies, pointers, q.size());

    // construct external forces & damping force
    externalForce.resize(q.size());
    externalForce.setZero();
    dampingForce.resize(qdot.size());
    dampingForce.setZero();
    // initialize damping force & damping matrix
    dampingForce = -damping * qdot;
    B.resize(qdot.size(), qdot.size());
    B.setIdentity();
    B *= -damping;

    // just add gravity to external forces
    Eigen::Vector3d g(0, -10, 0);
    Eigen::VectorXd gravity;
    gravity.resize(externalForce.size());
    gravity.setZero();
    for (int i = 0; i < q.size() / 3; i++) {
        gravity.segment(3*i, 3) += g;
    }
    externalForce = M * gravity;

    // initialize contact forces with boundaries
    boundaryForce.resize(externalForce.size());
    boundaryForce.setZero();

    // precompute M_inv * f_ext == total acceleration of external forces
    LUSolveM.analyzePattern(M);
    LUSolveM.factorize(M);
    Minv_fext.resize(externalForce.size());
    Minv_fext = LUSolveM.solve(externalForce); // in fact this is the same with gravity

    // set hash table
    double length = aveEdgeLength(softBodies);
    numX = int(6 / length);
    numY = int(6 / length);
    numZ = int(6 / length);
    gridSpacing = length;
    SpatialHash table(numX, numY, numZ, gridSpacing);

    // set boundaries
    Boundary floor; // default boundary is floor
    table.boundaryKeys(floor);
    boundaryList.push_back(floor);
    // Boundary leftWall(Eigen::Vector3d(-1, -1, -1), Eigen::Vector3d(1, 0, 0));
    // table.boundaryKeys(leftWall);
    // boundaryList.push_back(leftWall);

    // set boundary conditions
    std::vector<int> fixed;
    fixed.push_back(0);
    // fixed.push_back(1);
    // fixed.push_back(2);

    fixedPointIndexList.resize(fixed.size());
    // determine the exclusive constraint for a single object
    find_min_vertices(fixedPointIndexList[fixed[0]], softBodies[fixed[0]].V, 0.1);
    // find_min_vertices(fixedPointIndexList[fixed[1]], softBodies[fixed[1]].V, 0.1);

    // determine some objects with the same constraint
    // findMultiMinVertices(fixedPointIndexList, softBodies, fixed, 0.1);
    // findMultiMaxVertices(fixedPointIndexList, softBodies, fixed, 0.1);

    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixedMultiPointsConstraints(P, fixedPointIndexList, fixed, softBodies, pointers);

    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();

    std::cout << "finish regulate to the right size" << std::endl;

    // Visualize::viewer().callback_key_down = &keydown_callback;
}

int main(int argc, char **argv) {

    std::cout<<"Start\n";

    // assignment_setup(argc, argv, q, qdot);
    // initialize(q, qdot);
    initializeMultiBodySimulation(q, qdot);

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    //setup libigl viewer and activate
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();

    return 1;

}
