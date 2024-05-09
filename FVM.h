//Steady Finite Volume 2D Simulation, Rectangular mesh
#include<iostream>
#include "../../FEM/SIM-Tool/Core/eigen-3.4.0/Eigen/Dense"
#include<string>
#include<vector>
#include<cmath>
using namespace std;
using namespace Eigen;

//--------------------------------------DEFINE PARAMETERS AND METHOD----------------------------------------
#ifndef DEFINE_H
    #define DEFINE_H

    //Vertex type :
    #define TYPE_F 0
    #define TYPE_D 1
    #define TYPE_N 2
    #define TYPE_MIXED 3

    //Flow cases : 
    // #define Main_Problem_Task1
    // #define Versteeg_Example41  //Source free heat conduction in an insulated rod, validation for diffusion term            --Test for diffusion
    // #define Versteeg_Example42  //Heat conduction in an insulated rod, validation for diffusion term with source         --Test for source
    // #define Versteeg_Example51 //Convection diffusion flow                                                               --Test for convection
    #define Versteeg_FalseDiffusion //Convection flow for multidimensional flow

    //vertex properties : 
    static vector<double> vertex_coordinate_x;
    static vector<double> vertex_coordinate_y;
    static vector<short> vertex_BC; //BC type
    static vector<double> vertex_BC_value;  //BC value

    //dx and dy for vertex
    static vector<double> dx;   //dx of vertex (global)
    static vector<double> dy;   //dy of vertex (global)
    static vector<double> vertex_dx;    //dx of vertex (local)
    static vector<double> vertex_dy;    //dy of vertex (local)

    #if defined Main_Problem_Task1
        //Simulation Parameters : 
        static const double L1 = 1;    //Length y
        static const double L2 = 20;    //Length x
        static const int Nx = 5;        //Total number of vertex in x direction
        static const int Ny = 10;        //Total number of vertex in y direction
        static const int N_center = (Nx-1)*(Ny-1);
        static double phi_wb = 100;     //Assume uniform value of phi
        static double phi_sb = 50;      //Assume uniform value of phi, Dirichlet



        //dx and dy
        static const int num_dx = 1;
        static const int num_dy = 1;
        static double boundary_dx[num_dx] = {L2}; //written in x coordinate, rightmost coordinate
        static double boundary_dy[num_dy] = {L1}; //written in y coordinate, uppermost coordinate
        static int Nx_dx[num_dx] = {Nx};
        static int Ny_dy[num_dy] = {Ny};


        //Upper boundary division : 
        static const int num_BC_top = 1;
        static short type_top[num_BC_top] = {TYPE_N};     //Type of BC
        static double BC_coor_top[num_BC_top] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_top[num_BC_top] = {0};

        //Lower boundary division : 
        static const int num_BC_bottom = 1;
        static short type_bottom[num_BC_bottom] = {TYPE_D};     //Type of BC
        static double BC_coor_bottom[num_BC_bottom] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_bottom[num_BC_bottom] = {phi_sb};

        //Left boundary division : 
        static const int num_BC_left = 1;
        static short type_left[num_BC_left] = {TYPE_D};     //Type of BC
        static double BC_coor_left[num_BC_left] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_left[num_BC_top] = {phi_wb};

        //Right boundary division : 
        static const int num_BC_right = 1;
        static short type_right[num_BC_right] = {TYPE_N};    //Type of BC
        static double BC_coor_right[num_BC_right] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_right[num_BC_right] = {0};

        //Fluid properties : 
        static double rho = 1;
        static double vx = 25;
        static double vy = 5;
        static double Gamma_phi = 2;

        //Source term Q = a + b*phi
        static const double Source_Constant = 10;
        static const double Source_Linear = -2;
    #elif defined Versteeg_Example41
        //Simulation Parameters : 
        static const double L1 = 0.01;    //Length y
        static const double L2 = 0.5;    //Length x
        static const int Nx = 6;        //Total number of vertex in x direction
        static const int Ny = 2;        //Total number of vertex in y direction
        static const int N_center = (Nx-1)*(Ny-1);
        static double phi_wb = 100;     //Assume uniform value of phi
        static double phi_sb = 50;      //Assume uniform value of phi, Dirichlet

        



        //dx and dy
        static const int num_dx = 1;
        static const int num_dy = 1;
        static double boundary_dx[num_dx] = {L2}; //written in x coordinate, rightmost coordinate
        static double boundary_dy[num_dy] = {L1}; //written in y coordinate, uppermost coordinate
        static int Nx_dx[num_dx] = {Nx};
        static int Ny_dy[num_dy] = {Ny};
        


        //Upper boundary division : 
        static const int num_BC_top = 1;
        static short type_top[num_BC_top] = {TYPE_N};     //Type of BC
        static double BC_coor_top[num_BC_top] = {L2};     //BC boundary coordinate, defined over x coordinate
        static double BC_value_top[num_BC_top] = {0};

        //Lower boundary division : 
        static const int num_BC_bottom = 1;
        static short type_bottom[num_BC_bottom] = {TYPE_N};     //Type of BC
        static double BC_coor_bottom[num_BC_bottom] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_bottom[num_BC_bottom] = {0};

        //Left boundary division : 
        static const int num_BC_left = 1;
        static short type_left[num_BC_left] = {TYPE_D};     //Type of BC
        static double BC_coor_left[num_BC_left] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_left[num_BC_left] = {100};

        //Right boundary division : 
        static const int num_BC_right = 1;
        static short type_right[num_BC_right] = {TYPE_D};    //Type of BC
        static double BC_coor_right[num_BC_right] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_right[num_BC_right] = {500};

        //Fluid properties : 
        static double rho = 1;
        static double vx = 0;
        static double vy = 0;
        static double Gamma_phi = 1000;

        //Source term Q = a + b*phi
        static const double Source_Constant = 0;
        static const double Source_Linear = 0;
    #elif defined Versteeg_Example42
        //Simulation Parameters : 
        static const double L1 = 1;    //Length y
        static const double L2 = 0.02;    //Length x
        static const int Nx = 6;        //Total number of vertex in x direction
        static const int Ny = 2;        //Total number of vertex in y direction
        static const int N_center = (Nx-1)*(Ny-1);
        static double phi_wb = 100;     //Assume uniform value of phi
        static double phi_sb = 50;      //Assume uniform value of phi, Dirichlet

        



        //dx and dy
        static const int num_dx = 1;
        static const int num_dy = 1;
        static double boundary_dx[num_dx] = {L2}; //written in x coordinate, rightmost coordinate
        static double boundary_dy[num_dy] = {L1}; //written in y coordinate, uppermost coordinate
        static int Nx_dx[num_dx] = {Nx};
        static int Ny_dy[num_dy] = {Ny};
        


        //Upper boundary division : 
        static const int num_BC_top = 1;
        static short type_top[num_BC_top] = {TYPE_N};     //Type of BC
        static double BC_coor_top[num_BC_top] = {L2};     //BC boundary coordinate, defined over x coordinate
        static double BC_value_top[num_BC_top] = {0};

        //Lower boundary division : 
        static const int num_BC_bottom = 1;
        static short type_bottom[num_BC_bottom] = {TYPE_N};     //Type of BC
        static double BC_coor_bottom[num_BC_bottom] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_bottom[num_BC_bottom] = {0};

        //Left boundary division : 
        static const int num_BC_left = 1;
        static short type_left[num_BC_left] = {TYPE_D};     //Type of BC
        static double BC_coor_left[num_BC_left] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_left[num_BC_left] = {100};

        //Right boundary division : 
        static const int num_BC_right = 1;
        static short type_right[num_BC_right] = {TYPE_D};    //Type of BC
        static double BC_coor_right[num_BC_right] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_right[num_BC_right] = {200};

        //Fluid properties : 
        static double rho = 1;
        static double vx = 0;
        static double vy = 0;
        static double Gamma_phi = 0.5;

        //Source term Q = a + b*phi
        static const double Source_Constant = 1000*1000;
        static const double Source_Linear = 0;
    #elif defined Versteeg_Example51
        //Simulation Parameters : 
        static const double L1 = 1;    //Length y
        static const double L2 = 1;    //Length x
        static const int Nx = 6;        //Total number of vertex in x direction
        static const int Ny = 2;        //Total number of vertex in y direction
        static const int N_center = (Nx-1)*(Ny-1);
        static double phi_wb = 100;     //Assume uniform value of phi
        static double phi_sb = 50;      //Assume uniform value of phi, Dirichlet

        //dx and dy
        static const int num_dx = 1;
        static const int num_dy = 1;
        static double boundary_dx[num_dx] = {L2}; //written in x coordinate, rightmost coordinate
        static double boundary_dy[num_dy] = {L1}; //written in y coordinate, uppermost coordinate
        static int Nx_dx[num_dx] = {Nx};
        static int Ny_dy[num_dy] = {Ny};
        


        //Upper boundary division : 
        static const int num_BC_top = 1;
        static short type_top[num_BC_top] = {TYPE_N};     //Type of BC
        static double BC_coor_top[num_BC_top] = {L2};     //BC boundary coordinate, defined over x coordinate
        static double BC_value_top[num_BC_top] = {0};

        //Lower boundary division : 
        static const int num_BC_bottom = 1;
        static short type_bottom[num_BC_bottom] = {TYPE_N};     //Type of BC
        static double BC_coor_bottom[num_BC_bottom] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_bottom[num_BC_bottom] = {0};

        //Left boundary division : 
        static const int num_BC_left = 1;
        static short type_left[num_BC_left] = {TYPE_D};     //Type of BC
        static double BC_coor_left[num_BC_left] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_left[num_BC_left] = {1};

        //Right boundary division : 
        static const int num_BC_right = 1;
        static short type_right[num_BC_right] = {TYPE_D};    //Type of BC
        static double BC_coor_right[num_BC_right] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_right[num_BC_right] = {0};

        //Fluid properties : 
        static double rho = 1;
        static double vx = 2.5;
        static double vy = 0;
        static double Gamma_phi = 0.1;

        //Source term Q = a + b*phi
        static const double Source_Constant = 0;
        static const double Source_Linear = 0;
    #elif defined Versteeg_FalseDiffusion
        //Simulation Parameters : 
        static const double L1 = 1;    //Length y
        static const double L2 = 1;    //Length x
        static const int Nx = 11;        //Total number of vertex in x direction
        static const int Ny = 11;        //Total number of vertex in y direction
        static const int N_center = (Nx-1)*(Ny-1);
        static double phi_wb = 100;     //Assume uniform value of phi
        static double phi_sb = 50;      //Assume uniform value of phi, Dirichlet

        //dx and dy
        static const int num_dx = 1;
        static const int num_dy = 1;
        static double boundary_dx[num_dx] = {L2}; //written in x coordinate, rightmost coordinate
        static double boundary_dy[num_dy] = {L1}; //written in y coordinate, uppermost coordinate
        static int Nx_dx[num_dx] = {Nx};
        static int Ny_dy[num_dy] = {Ny};
        


        //Upper boundary division : 
        static const int num_BC_top = 1;
        static short type_top[num_BC_top] = {TYPE_D};     //Type of BC
        static double BC_coor_top[num_BC_top] = {L2};     //BC boundary coordinate, defined over x coordinate
        static double BC_value_top[num_BC_top] = {100};

        //Lower boundary division : 
        static const int num_BC_bottom = 1;
        static short type_bottom[num_BC_bottom] = {TYPE_D};     //Type of BC
        static double BC_coor_bottom[num_BC_bottom] = {L2};     //BC coordinate, defined over x coordinate
        static double BC_value_bottom[num_BC_bottom] = {0};

        //Left boundary division : 
        static const int num_BC_left = 1;
        static short type_left[num_BC_left] = {TYPE_D};     //Type of BC
        static double BC_coor_left[num_BC_left] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_left[num_BC_left] = {100};

        //Right boundary division : 
        static const int num_BC_right = 1;
        static short type_right[num_BC_right] = {TYPE_D};    //Type of BC
        static double BC_coor_right[num_BC_right] = {L1};     //BC coordinate, defined over y coordinate
        static double BC_value_right[num_BC_right] = {0};

        //Fluid properties : 
        static double rho = 1;
        static double vx = 2;
        static double vy = 2;
        static double Gamma_phi = 0;

        //Source term Q = a + b*phi
        static const double Source_Constant = 0;
        static const double Source_Linear = 0;
    #endif
    
    //Diffusive flux method :
    #define Diffusive_CentralDifference

    //Convective flux method : 
    #define Convective_CentralDifference
    // #define Convective_UpwindScheme
    // #define Convective_SOU
    // #define Convective_QUICK
    // #define Convective_FROMM

    //Method for finding neighbor : 
    #define General_Method  
    //The general method states that any cell center is considered as neighbor if they share the same face
    //2 point vertex in 2d (line), 4 point vertex in 3d (face)

#endif

//-----------------------------------------------------------------------------------------------------------




//-------------------------------------------CLASS FOR SIMULATION--------------------------------------------
#ifndef CLASS_H
    class CenterPoint {
        public : 
            //Variables :
            //Center coordinate 
            double center_coordinate[2]; //(x,y)
            int center_index[2]; //(i,j)
            int position_global; //Number index of center

            //Neighbor number
            vector<int> neighbor_number;    //E, N, W, S CCW direction, in number not coordinate

            #ifdef General_Method
                vector<Vector2i> shared_face;
                vector<Vector3d> Normal;
            #endif

            //For nonuniform grid :
            vector<double> dx;      //E, N, W, S CCW direction
            vector<double> dy;      //E, N, W, S CCW direction

            //Vertex number
            vector<int> vertex_number;   //Given in number not coordinate
            vector<double> vertex_coordinate_x; 
            vector<double> vertex_coordinate_y; 
            vector<short> vertex_type;
            vector<double> vertex_BC_value;

            //Fluid
            double density;
            double velocity[2];
            double phi;
            double Gamma_phi;
            double Peclet_x;
            double Peclet_y;

            //Finite volume property : 
            double a_east = 0;
            double a_west = 0;
            double a_north = 0;
            double a_south = 0;
            double a_center = 0;
            

            //Mutator : 
            void Find_Center();     //Given Vertex, find the center of the cell 
            void Find_Neighbor();   //Given the center of the cell, find the neighbor
    };

    class FVM {
        public : 
            //Variables : 
            CenterPoint*center;
            MatrixXd A = MatrixXd::Zero(N_center, N_center);
            MatrixXd b = MatrixXd::Zero(N_center, 1);
            MatrixXd phi = MatrixXd::Zero(N_center, 1);

            //Mutator :
            void Initialization();
            void Local2Global();
            void Finite_Volume_Method();    //Finding a's coefficient, treatment of boundary condition and finding phi
    };

#endif

//-----------------------------------------------------------------------------------------------------------





//-----------------------------------------------OUTPUT FILES------------------------------------------------
#ifndef OUTPUT_H
    #define OUTPUT_H

    void PrintInfo();
    void PrintCenter(FVM& fvm);
    #ifdef General_Method
        int Linear_Search(vector<int> target, Vector2i Wanted);
        void Find_Neighbor(FVM& fvm);
    #endif

    //Output CSV data : 
    void OutputCSV(FVM& fvm, string name);
#endif
//-----------------------------------------------------------------------------------------------------------