//Steady Finite Volume 2D Simulation, Rectangular mesh
#include<iostream>
#include "../../FEM/SIM-Tool/Core/eigen-3.4.0/Eigen/Dense"
#include <string>
#include "FVM.h"
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;
using namespace Eigen;

//Definition of member function in CenterPoint class :
void CenterPoint::Find_Center() {   
    int i_center = position_global%(Nx-1), j_center = position_global/(Nx-1);
    int i_vertex = i_center + 0.5, j_vertex = j_center + 0.5;
    center_index[0] = i_center; center_index[1] = j_center;

    double sum_x = 0, sum_y = 0;
    for(int i = 0; i<vertex_number.size(); i++) {
        sum_x += vertex_coordinate_x[i];
        sum_y += vertex_coordinate_y[i];
    }
    center_coordinate[0] = sum_x/vertex_number.size(); center_coordinate[1] = sum_y/vertex_number.size();
    
}

void CenterPoint::Find_Neighbor() {
    //Given the center of the cell, find the neighbor (for structured grid)
    #if defined General_Method  //For general grid
        //Finding the set of shared face and finding normal :
        Vector3d vec1(0,0,1);
        for(int n = 0; n<vertex_number.size(); n++) {
            int idx1 = n%vertex_number.size(), idx2 = (n+1)%vertex_number.size();
            Vector2i shared_faces(vertex_number[idx1], vertex_number[idx2]);
            shared_face.push_back(shared_faces);

            int vertex1_index = vertex_number[idx1], vertex2_index = vertex_number[idx2];
            double x1 = vertex_coordinate_x[n], y1 = vertex_coordinate_y[n];
            double x2 = vertex_coordinate_x[(n+1)%4], y2 = vertex_coordinate_y[(n+1)%4];

            Vector3d vec2(x2-x1, y2-y1, 0);
            Vector3d normal = vec2.cross(vec1);
            double norm = sqrt(normal.dot(normal));
            Normal.push_back(normal/norm);
        }

        
        

    #endif
}

//Definition of member function in FVM class : 
void FVM::Initialization() {
    //Memory allocation for center:
    center = new CenterPoint [N_center];

    //Vertex part : 
    //1. dx and dy for vertex
    //a. dx : 
    for(int i = 0; i<num_dx; i++) {
        double x1, x2;
        if(i == 0) {
            x1 = 0; x2 = boundary_dx[i];
        }
        else {
            x1 = boundary_dx[i-1]; x2 = boundary_dx[i];
        }
        dx.push_back((x2-x1)/(Nx_dx[i]-1)); //For global purpose

    }

    //b. dy : 
    for(int i = 0; i<num_dy; i++) {
        double y1, y2;
        if(i == 0) {
            y1 = 0; y2 = boundary_dy[i];
        }
        else {
            y1 = boundary_dy[i-1]; y2 = boundary_dy[i];
        }
        dy.push_back((y2-y1)/(Ny_dy[i]-1));   //For global purpose
    }
    
    //2. Coordinate and BC for vertex : 
    for(int vertex_num = 0; vertex_num<Nx*Ny; vertex_num++) {
        //a. x coordinate : 
        double x = 0;
        int idx_i = vertex_num%Nx;
        int idx_dx_trial = Nx_dx[0]-1;
        for(int idx_dx = 0; idx_dx<num_dx; idx_dx++) {
            if(idx_i > idx_dx_trial) {
                x += (Nx_dx[idx_dx] - 1)*dx[idx_dx];
                idx_dx_trial += Nx_dx[idx_dx+1];
            }
            else {
                if(idx_dx == 0) {
                    x += (idx_i)*dx[idx_dx];
                }
                else {
                    x += (idx_i - Nx_dx[idx_dx-1] + 1 )*dx[idx_dx];
                }
                vertex_dx.push_back(dx[idx_dx]);
                vertex_coordinate_x.push_back(x);                
                idx_dx = num_dx+1;
                
            }
        }

        //b. y coordinate : 
        double y = 0;
        int idx_j = int(vertex_num/Nx);
        int idx_dy_trial = Ny_dy[0]-1;
        for(int idx_dy = 0; idx_dy<num_dy; idx_dy++) {
            if(idx_j > idx_dy_trial) {
                y += (Ny_dy[idx_dy] - 1)*dy[idx_dy];
                idx_dy_trial += Ny_dy[idx_dy+1];
            }
            else {
                if(idx_dy == 0) {
                    y += (idx_j)*dy[idx_dy];
                }
                else {
                    y += (idx_j - Ny_dy[idx_dy-1] + 1 )*dy[idx_dy];
                }
                vertex_dy.push_back(dy[idx_dy]);
                vertex_coordinate_y.push_back(y);
                idx_dy = num_dy+1;
            }
        }

        //Boundary condition : 
        if(x == 0) {
            //Consider left BC : 
            for(int j = 0; j<num_BC_left; j++) {
                if(y <= BC_coor_left[j]) {
                    vertex_BC.push_back(type_left[j]);
                    vertex_BC_value.push_back(BC_coor_left[j]);
                    j = num_BC_left;
                }
            }
        }
        else if (x == L2) {
            //Consider right BC : 
            for(int j = 0; j<num_BC_right; j++) {
                if(y <= BC_coor_right[j]) {
                    vertex_BC.push_back(type_right[j]);
                    vertex_BC_value.push_back(BC_coor_right[j]);
                    j = num_BC_right;
                }
            }
        }
        else if (y == 0) {
            //Consider bottom BC : 
            for(int i = 0; i<num_BC_bottom; i++) {
                if(x <= BC_coor_bottom[i]) {
                    vertex_BC.push_back(type_bottom[i]);
                    vertex_BC_value.push_back(BC_coor_bottom[i]);
                    i = num_BC_bottom;
                }
            }
        }
        else if (y == L1) {
            //Consider top BC : 
            for(int i = 0; i<num_BC_top; i++) {
                if(x <= BC_coor_top[i]) {
                    vertex_BC.push_back(type_top[i]);
                    vertex_BC_value.push_back(BC_coor_top[i]);
                    i = num_BC_top;
                }
            }
        }
        else {
            vertex_BC.push_back(TYPE_F);
        }
    }


    //Last : Inputing vertex coordinate to center
    for(int num = 0; num<N_center; num++) {
        int position_global = num; 
        int i_center1 = num%(Nx-1), j_center1 = num/(Nx-1);   //in center coordinate : 
        double i_center2 = i_center1 + 0.5, j_center2 = j_center1 + 0.5;    //In vertex coordinate : 

        int i_vertex, j_vertex, num_vertex;
        double x, y, dx, dy;
        
        //Vertex 0 : 
        i_vertex = ceil(i_center2); j_vertex = floor(j_center2);
        num_vertex = j_vertex*Nx+i_vertex;
        x = vertex_coordinate_x[num_vertex]; y = vertex_coordinate_y[num_vertex]; dx = vertex_dx[num_vertex]; dy = vertex_dy[num_vertex];
        center[num].vertex_number.push_back(num_vertex);
        center[num].vertex_coordinate_x.push_back(x);
        center[num].vertex_coordinate_y.push_back(y);
        center[num].vertex_type.push_back(vertex_BC[num_vertex]);
        center[num].vertex_BC_value.push_back(vertex_BC_value[num_vertex]);

        //Vertex 1 : 
        i_vertex = ceil(i_center2); j_vertex = ceil(j_center2);
        num_vertex = j_vertex*Nx+i_vertex;
        x = vertex_coordinate_x[num_vertex]; y = vertex_coordinate_y[num_vertex]; dx = vertex_dx[num_vertex]; dy = vertex_dy[num_vertex];
        center[num].vertex_number.push_back(num_vertex);
        center[num].vertex_coordinate_x.push_back(x);
        center[num].vertex_coordinate_y.push_back(y);
        center[num].vertex_type.push_back(vertex_BC[num_vertex]);
        center[num].vertex_BC_value.push_back(vertex_BC_value[num_vertex]);

        //Vertex 2 : 
        i_vertex = floor(i_center2); j_vertex = ceil(j_center2);
        num_vertex = j_vertex*Nx+i_vertex;
        x = vertex_coordinate_x[num_vertex]; y = vertex_coordinate_y[num_vertex]; dx = vertex_dx[num_vertex]; dy = vertex_dy[num_vertex];
        center[num].vertex_number.push_back(num_vertex);
        center[num].vertex_coordinate_x.push_back(x);
        center[num].vertex_coordinate_y.push_back(y);
        center[num].vertex_type.push_back(vertex_BC[num_vertex]);
        center[num].vertex_BC_value.push_back(vertex_BC_value[num_vertex]);

        //Vertex 3 : 
        i_vertex = floor(i_center2); j_vertex = floor(j_center2);
        num_vertex = j_vertex*Nx+i_vertex;
        x = vertex_coordinate_x[num_vertex]; y = vertex_coordinate_y[num_vertex]; dx = vertex_dx[num_vertex]; dy = vertex_dy[num_vertex];
        center[num].vertex_number.push_back(num_vertex);
        center[num].vertex_coordinate_x.push_back(x);
        center[num].vertex_coordinate_y.push_back(y);
        center[num].vertex_type.push_back(vertex_BC[num_vertex]);
        center[num].vertex_BC_value.push_back(vertex_BC_value[num_vertex]);
    }

    //Finding center and neighbor
    for(int i = 0; i<N_center; i++) {
        center[i].position_global = i;
        center[i].Find_Center();
        center[i].Find_Neighbor();
    }

    //Material properties and fluid properties : 
    for(int i = 0; i<N_center; i++) {
        //Vary the fluid properties here
        center[i].density = rho;
        center[i].velocity[0] = vx;
        center[i].velocity[1] = vy;
        center[i].Gamma_phi = Gamma_phi;
    }

}



void FVM::Local2Global() {  
    for(int i = 0; i<N_center; i++) {
        //FVM coefficients : 
        double a[5] = {center[i].a_center, center[i].a_east, center[i].a_north, center[i].a_west, center[i].a_south};
        int idx[5] = {i,center[i].neighbor_number[0], center[i].neighbor_number[1], center[i].neighbor_number[2], center[i].neighbor_number[3]};

        
        for(int j = 0; j<5; j++) {
            if(idx[j] != -1) {A(i,idx[j]) = a[j];}
        }
    }
}

void FVM::Finite_Volume_Method() {
    for(int i = 0; i<N_center; i++) {
        //Calculate east, west, north and south properties (for nonuniform grid) : 
        double g_ratio[center[i].vertex_number.size()]; //e, n, w, s counterclockwise direction
        double Gamma[center[i].vertex_number.size()];   //e, n, w, s counterclockwise direction
        double rho[center[i].vertex_number.size()];
        double velocity_x[center[i].vertex_number.size()];
        double velocity_y[center[i].vertex_number.size()];
        double distance_Cface, distance_Cneighbor;

        for(int n = 0; n<center[i].vertex_number.size(); n++) {
            switch(n) {
                case 0 :
                case 2 :distance_Cface = abs(vertex_coordinate_x[center[i].vertex_number[n]] - center[i].center_coordinate[0]);
                        if(center[i].neighbor_number[n]!=-1) {
                            distance_Cneighbor = abs(center[i].center_coordinate[0] - center[center[i].neighbor_number[n]].center_coordinate[0]);
                            g_ratio[n] = distance_Cface/distance_Cneighbor; 
                        }
                        else{g_ratio[n] = 0.5;}
                        
                        break;
                case 1 : 
                case 3 :distance_Cface = abs(vertex_coordinate_y[center[i].vertex_number[n]] - center[i].center_coordinate[1]);
                        if(center[i].neighbor_number[n]!=-1) {
                            distance_Cneighbor = abs(center[i].center_coordinate[1] - center[center[i].neighbor_number[n]].center_coordinate[1]);
                            g_ratio[n] = distance_Cface/distance_Cneighbor;
                        }
                        else {g_ratio[n] = 0.5;}
                        
                        break;
            }
            
            double Gamma_neighbor, rho_neighbor, velocityx_neighbor, velocityy_neighbor;
            if(center[i].neighbor_number[n] != -1) {
                Gamma_neighbor = center[i].Gamma_phi;
                rho_neighbor = center[i].density;
                velocityx_neighbor = center[i].velocity[0];
                velocityy_neighbor = center[i].velocity[1];
            }
            else{
                Gamma_neighbor = center[n].Gamma_phi;
                rho_neighbor = center[n].density;
                velocityx_neighbor = center[n].velocity[0];
                velocityy_neighbor = center[n].velocity[1];
            }
            
            double Gamma_center = center[n].Gamma_phi;
            double rho_center = center[n].density;
            double velocityx_center = center[n].velocity[0];
            double velocityy_center = center[n].velocity[1];
            Gamma[n] = (1-g_ratio[n])*Gamma_center + g_ratio[n]*Gamma_neighbor; 
            rho[n] = (1-g_ratio[n])*rho_center + g_ratio[n]*rho_neighbor; 
            velocity_x[n] = (1-g_ratio[n])*velocityx_center + g_ratio[n]*velocityx_neighbor; 
            velocity_y[n] = (1-g_ratio[n])*velocityy_center + g_ratio[n]*velocityy_neighbor; 
        }
        
        
        //Calculate a_east, a_west, a_north, a_south, a_center and b:
        //Diffusion term : 
        #if defined Diffusive_CentralDifference and defined Convective_UpwindScheme
            for(int n = 0; n<center[i].neighbor_number.size(); n++) {
                Vector3d velocity_neighbor(velocity_x[n], velocity_y[n], 0);
                double v_dot_n = velocity_neighbor.dot(center[i].   Normal[n]);
                double mdot_neighbor;
                
                //Calculate mdot : 
                switch(n) {
                    case 0 : 
                    case 2 : mdot_neighbor = rho[n]*v_dot_n*center[i].dy[n]; break;
                    case 1 : 
                    case 3 : mdot_neighbor = rho[n]*v_dot_n*center[i].dx[n]; break;
                }

                //Calculate the center coefficients : 
                if(n%2 == 0) {
                    center[i].a_center += -(Gamma[n]*center[i].dy[n]/center[i].dx[n] + max(mdot_neighbor,0.));
                }
                else {
                    center[i].a_center += -(Gamma[n]*center[i].dx[n]/center[i].dy[n] + max(mdot_neighbor,0.));
                }

                //Calculate the face coefficients : 
                int neighbor_number = center[i].neighbor_number[n];
                if(neighbor_number != -1) {
                    switch(n) {
                        case 0 : center[i].a_east   = Gamma[n]*center[i].dy[n]/center[i].dx[n] + max(-mdot_neighbor,0.); break;
                        case 1 : center[i].a_north  = Gamma[n]*center[i].dx[n]/center[i].dy[n] + max(-mdot_neighbor,0.); break;
                        case 2 : center[i].a_west   = Gamma[n]*center[i].dy[n]/center[i].dx[n] + max(-mdot_neighbor,0.); break;
                        case 3 : center[i].a_south  = Gamma[n]*center[i].dx[n]/center[i].dy[n] + max(-mdot_neighbor,0.); break;
                    }
                }
                else {
                    //Check for boundary condition : 
                    double phi_boundary, q;
                    short boundary_type;
                    if(n == 0) {
                        double y_BC = center[i].center_coordinate[1];
                        //Check for right BC : 
                        for(int idxBC = 0; idxBC<num_BC_right; idxBC++) {
                            if(y_BC<=BC_coor_right[idxBC]) {
                                phi_boundary = BC_value_right[idxBC];
                                q = BC_value_right[idxBC];
                                boundary_type = type_right[idxBC];
                                idxBC = num_BC_right;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_east = Gamma[n]*center[i].dy[n]/center[i].dx[n] + max(-mdot_neighbor,0.);
                        if(boundary_type == TYPE_D) {
                            center[i].a_east = 0;
                            b(i) += -a_east*phi_boundary/g_ratio[n];
                            center[i].a_center += a_east*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_east = 0;
                            b(i) += q*center[i].dx[n]/Gamma[n]*a_east;
                            center[i].a_center += a_east;
                        }
                    }
                    else if (n == 1) {
                        double x_BC = center[i].center_coordinate[0];
                        //Check for top BC : 
                        for(int idxBC = 0; idxBC<num_BC_top; idxBC++) {
                            if(x_BC<=BC_coor_top[idxBC]) {
                                phi_boundary = BC_value_top[idxBC];
                                q = BC_value_top[idxBC];
                                boundary_type = type_top[idxBC];
                                idxBC = num_BC_top;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_north = Gamma[n]*center[i].dx[n]/center[i].dy[n] + max(-mdot_neighbor,0.);
                        if(boundary_type == TYPE_D) {
                            center[i].a_north = 0;
                            b(i) += -a_north*phi_boundary/g_ratio[n];
                            center[i].a_center += a_north*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_north = 0;
                            b(i) += q*center[i].dy[n]/Gamma[n]*(-Gamma[n]*center[i].dx[n]/center[i].dy[n]);
                            center[i].a_center += a_north;
                        }
                    }
                    else if (n == 2) {
                        double y_BC = center[i].center_coordinate[1];
                        //Check for left BC : 
                        for(int idxBC = 0; idxBC<num_BC_left; idxBC++) {
                            if(y_BC<=BC_coor_left[idxBC]) {
                                phi_boundary = BC_value_left[idxBC];
                                q = BC_value_left[idxBC];
                                boundary_type = type_left[idxBC];
                                idxBC = num_BC_right;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_west = Gamma[n]*center[i].dy[n]/center[i].dx[n] + max(-mdot_neighbor,0.);
                        if(boundary_type == TYPE_D) {
                            center[i].a_west = 0;
                            b(i) += -a_west*phi_boundary/g_ratio[n];
                            center[i].a_center += a_west*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_west = 0;
                            b(i) += q*center[i].dx[n]/Gamma[n]*(-Gamma[n]*center[i].dy[n]/center[i].dx[n]);
                            center[i].a_center += a_west;
                        }
                    }
                    else if (n == 3) {
                        double x_BC = center[i].center_coordinate[0];
                        //Check for bottom BC : 
                        for(int idxBC = 0; idxBC<num_BC_bottom; idxBC++) {
                            if(x_BC<=BC_coor_bottom[idxBC]) {
                                phi_boundary = BC_value_bottom[idxBC];
                                q = BC_value_bottom[idxBC];
                                boundary_type = type_bottom[idxBC];
                                idxBC = num_BC_bottom;
                            }
                        }

                        //look at the type of BC : 
                        double a_south = Gamma[n]*center[i].dx[n]/center[i].dy[n] + max(-mdot_neighbor,0.);
                        if(boundary_type == TYPE_D) {
                            center[i].a_south = 0;
                            b(i) += -a_south*phi_boundary/g_ratio[n];
                            center[i].a_center += a_south*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_south = 0;
                            b(i) += q*center[i].dy[n]/Gamma[n]*(-Gamma[n]*center[i].dx[n]/center[i].dy[n]);
                            center[i].a_center += a_south;
                        }   
                    }
                }
                
            }
            
        #elif defined Diffusive_CentralDifference and defined Convective_CentralDifference
            for(int n = 0; n<center[i].neighbor_number.size(); n++) {
                Vector3d velocity_neighbor(velocity_x[n], velocity_y[n], 0);
                double v_dot_n = velocity_neighbor.dot(center[i].   Normal[n]);
                double mdot_neighbor;
                
                //Calculate mdot : 
                switch(n) {
                    case 0 : mdot_neighbor = rho[n]*v_dot_n*center[i].dy[n]; break;
                    case 1 : mdot_neighbor = rho[n]*v_dot_n*center[i].dx[n]; break;
                    case 2 : mdot_neighbor = rho[n]*v_dot_n*center[i].dy[n]; break;
                    case 3 : mdot_neighbor = rho[n]*v_dot_n*center[i].dx[n]; break;
                }

                //Calculate the center coefficients : 
                switch(n) {
                    case 0 : center[i].a_center += (Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2); break;
                    case 1 : center[i].a_center += (Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2); break;
                    case 2 : center[i].a_center += (Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2); break;
                    case 3 : center[i].a_center += (Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2); break;
                }


                //Calculate the face coefficients : 
                int neighbor_number = center[i].neighbor_number[n];
                if(neighbor_number != -1) {
                    switch(n) {
                        case 0 : center[i].a_east   = - Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2; break;
                        case 1 : center[i].a_north  = - Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2; break;
                        case 2 : center[i].a_west   = - Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2; break;
                        case 3 : center[i].a_south  = - Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2; break;
                    }
                }
                else {
                    //Check for boundary condition : 
                    double phi_boundary, q;
                    short boundary_type;
                    if(n == 0) {
                        double y_BC = center[i].center_coordinate[1];
                        //Check for right BC : 
                        for(int idxBC = 0; idxBC<num_BC_right; idxBC++) {
                            if(y_BC<=BC_coor_right[idxBC]) {
                                phi_boundary = BC_value_right[idxBC];
                                q = BC_value_right[idxBC];
                                boundary_type = type_right[idxBC];
                                idxBC = num_BC_right;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_east = - Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2;
                        if(boundary_type == TYPE_D) {
                            center[i].a_east = 0;
                            b(i) += -a_east*phi_boundary/g_ratio[n];
                            center[i].a_center += a_east*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_east = 0;
                            b(i) += q*center[i].dx[n]/Gamma[n]*a_east;
                            center[i].a_center += a_east;
                        }
                    }
                    else if (n == 1) {
                        double x_BC = center[i].center_coordinate[0];
                        //Check for top BC : 
                        for(int idxBC = 0; idxBC<num_BC_top; idxBC++) {
                            if(x_BC<=BC_coor_top[idxBC]) {
                                phi_boundary = BC_value_top[idxBC];
                                q = BC_value_top[idxBC];
                                boundary_type = type_top[idxBC];
                                idxBC = num_BC_top;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_north = - Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2;
                        if(boundary_type == TYPE_D) {
                            center[i].a_north = 0;
                            b(i) += -a_north*phi_boundary/g_ratio[n];
                            center[i].a_center += a_north*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_north = 0;
                            b(i) += q*center[i].dy[n]/Gamma[n]*(-Gamma[n]*center[i].dx[n]/center[i].dy[n]);
                            center[i].a_center += a_north;
                        }
                    }
                    else if (n == 2) {
                        double y_BC = center[i].center_coordinate[1];
                        //Check for left BC : 
                        for(int idxBC = 0; idxBC<num_BC_left; idxBC++) {
                            if(y_BC<=BC_coor_left[idxBC]) {
                                phi_boundary = BC_value_left[idxBC];
                                q = BC_value_left[idxBC];
                                boundary_type = type_left[idxBC];
                                idxBC = num_BC_right;
                            }
                        }
                        
                        //look at the type of BC : 
                        double a_west = - Gamma[n]*center[i].dy[n]/center[i].dx[n] + mdot_neighbor/2;
                        if(boundary_type == TYPE_D) {
                            center[i].a_west = 0;
                            b(i) += -a_west*phi_boundary/g_ratio[n];
                            center[i].a_center += a_west*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_west = 0;
                            b(i) += q*center[i].dx[n]/Gamma[n]*(-Gamma[n]*center[i].dy[n]/center[i].dx[n]);
                            center[i].a_center += a_west;
                        }
                    }
                    else if (n == 3) {
                        double x_BC = center[i].center_coordinate[0];
                        //Check for bottom BC : 
                        for(int idxBC = 0; idxBC<num_BC_bottom; idxBC++) {
                            if(x_BC<=BC_coor_bottom[idxBC]) {
                                phi_boundary = BC_value_bottom[idxBC];
                                q = BC_value_bottom[idxBC];
                                boundary_type = type_bottom[idxBC];
                                idxBC = num_BC_bottom;
                            }
                        }

                        //look at the type of BC : 
                        double a_south = - Gamma[n]*center[i].dx[n]/center[i].dy[n] + mdot_neighbor/2;
                        if(boundary_type == TYPE_D) {
                            center[i].a_south = 0;
                            b(i) += -a_south*phi_boundary/g_ratio[n];
                            center[i].a_center += a_south*(1-1/g_ratio[n]);
                        }
                        else if (boundary_type == TYPE_N) {
                            center[i].a_south = 0;
                            b(i) += q*center[i].dy[n]/Gamma[n]*(-Gamma[n]*center[i].dx[n]/center[i].dy[n]);
                            center[i].a_center += a_south;
                        }   
                    }
                }
                
            }
            

        #elif defined Diffusive_CentralDifference and defined Convective_SOU
        #elif defined Diffusive_CentralDifference and defined Convective_QUICK
        #elif defined Diffusive_CentralDifference and defined Convective_FROMM
        #endif

        //Source term : 
        //Center coefficient :
        double dx_vertex = abs(vertex_coordinate_x[center[i].vertex_number[3]] - vertex_coordinate_x[center[i].vertex_number[0]]), dy_vertex = abs(vertex_coordinate_y[center[i].vertex_number[1]] - vertex_coordinate_y[center[i].vertex_number[0]]);
        center[i].a_center += Source_Linear*dx_vertex*dy_vertex;

        //b coefficient : 
        b(i) += -Source_Constant*dx_vertex*dy_vertex;
    }

    //From local to global matrix : 
    Local2Global();

    //Matrix calculation to phi
    cout<<"Matrix A = \n"<<A<<endl;
    cout<<"Matrix b = \n"<<b<<endl;
    phi = A.colPivHouseholderQr().solve(b);
    cout<<"\nSolution phi = \n"<<phi<<endl;

    //Moving data to center : 
    for(int i = 0; i<N_center; i++) {
        center[i].phi = phi(i);
    }
}



//Other function definition : 
void PrintInfo() {
    cout<<"Case                                                    : 2D transport through rectangular duct with source generation\n";
    cout<<"Length x                                                : "<<L2<<endl;
    cout<<"Length y                                                : "<<L1<<endl;
    cout<<"Number of grids x                                       : "<<Nx<<endl;
    cout<<"Number of grids y                                       : "<<Ny<<"\n\n";

    cout<<"\nBoundary condition                         : \n";
    cout<<"       - Left                                           : Dirichlet (phi = "<<phi_wb<<")"<<endl;
    cout<<"       - Bottom                                         : Dirichlet (phi = "<<phi_sb<<")"<<endl;
    cout<<"       - Top                                            : Neumann (dphi/dy = 0)"<<endl;
    cout<<"       - Right                                          : Neumann (dphi/dx = 0)"<<endl;

    cout<<"\n\n";
}

void PrintCenter(FVM& fvm) {
    for(int i = 0; i<N_center; i++) {
        cout<<"Center - "<<fvm.center[i].position_global<<endl;
        cout<<"    center index      = ("<<fvm.center[i].center_index[0]<<", "<<fvm.center[i].center_index[1]<<")\n";
        cout<<"    center coordinate = ("<<fvm.center[i].center_coordinate[0]<<", "<<fvm.center[i].center_coordinate[1]<<")\n\n";
        
        cout<<"Vertex : \n";
        for(int num = 0; num<fvm.center[i].vertex_number.size(); num++) {
            cout<<fvm.center[i].vertex_number[num]<<", ";
        }

        cout<<"\nSet of (dx,dy) : \n";
        for(int num = 0; num<fvm.center[i].vertex_number.size(); num++) {
            cout<<"("<<fvm.center[i].dx[num]<<", "<<fvm.center[i].dy[num]<<"), ";
        }

        cout<<"\nvertex BC : \n";
        for(int num = 0; num<fvm.center[i].vertex_number.size(); num++) {
            if(fvm.center[i].vertex_type[num] == 0) {
                cout<<"Fluid, ";
            }
            else if(fvm.center[i].vertex_type[num] == 1) {
                cout<<"Dirichlet, ";
            }
            else if(fvm.center[i].vertex_type[num] == 2) {
                cout<<"Neumann, ";
            }
            else if(fvm.center[i].vertex_type[num] == 3) {
                cout<<"Mixed, ";
            }
        }

        // #ifdef General_Method
        //     cout<<"\nNormal : \n";
        //     for(int num = 0; num<fvm.center[i].Normal.size(); num++) {
        //         cout<<fvm.center[i].Normal[num]<<"\n\n";
        //     }
        // #endif

        cout<<"\nCenter Neighbor : \n";
        for(int num = 0; num<fvm.center[i].neighbor_number.size(); num++) {
            cout<<fvm.center[i].neighbor_number[num]<<", ";
        }
        cout<<"\n\n";

    }
}

void OutputCSV(FVM& fvm, string name) { 
	ofstream ofs;
	ofs.open(name);
	ofs << "x,y,density,velocity_x,velocity_y,phi\n";
	for (int i = 0; i<N_center; i++) {
        ofs <<fvm.center[i].center_coordinate[0]<< "," <<fvm.center[i].center_coordinate[1]<< ","<<fvm.center[i].density<< ","<<fvm.center[i].velocity[0]<< ","<<fvm.center[i].velocity[1]<< ","<<fvm.center[i].phi<<"\n";		
	}
	ofs.close();
}


#ifdef General_Method
    int Linear_Search(vector<int> search_target, Vector2i Wanted) {
        int matches = 0;
        // cout<<"Wanted = "<<Wanted[0]<<", "<<Wanted[1]<<endl;
        // cout<<"Search target = "<<search_target[0]<<", "<<search_target[1]<<", "<<search_target[2]<<", "<<search_target[3]<<"\n";
        for(int want_idx = 0; want_idx<2; want_idx++) {
            for(int idx = 0; idx<search_target.size(); idx++) {

                if(Wanted[want_idx] == search_target[idx]) {
                    matches++; idx = search_target.size()+1;
                }
            }
        }
        // cout<<"Matches = "<<matches<<"\n\n";
        return matches;
    }

    void Find_Neighbor(FVM& fvm) {
        for(int i = 0; i<N_center; i++) {
            for(int num = 0; num<fvm.center[i].shared_face.size(); num++) {
                Vector2i shared = fvm.center[i].shared_face[num];

                //Compare the wanted shared face to other and take the index
                int number_vertex = 2;
                int index_j = -1;
                for(int j = 0; j<N_center; j++) {
                    if(j!=i) {
                        // cout<<"j = "<<j<<", i = "<<i<<endl;
                        vector<int> set_vertex = fvm.center[j].vertex_number;

                        //Perform linear search on each vertex number
                        int match = Linear_Search(set_vertex, shared);
                        if(match == number_vertex) index_j = j;
                    }
                }

                fvm.center[i].neighbor_number.push_back(index_j);
            }

            //Finding dx and dy for nonuniform grid : 
            int neighbor_number;
            int x, y;
            for(int n = 0; n<fvm.center[i].vertex_number.size(); n++) {
                neighbor_number = fvm.center[i].neighbor_number[n];
                
                if(neighbor_number == -1) {
                    double x1 = fvm.center[i].vertex_coordinate_x[0], x2 = fvm.center[i].vertex_coordinate_x[3];
                    double y1 = fvm.center[i].vertex_coordinate_y[0], y2 = fvm.center[i].vertex_coordinate_y[1];
                    fvm.center[i].dx.push_back(abs(x2-x1)); 
                    fvm.center[i].dy.push_back(abs(y2-y1));
                }
                else {
                    if(n%2 == 0) {
                        //For x : 
                        double x_center = fvm.center[i].center_coordinate[0];
                        double x = fvm.center[neighbor_number].center_coordinate[0]; 
                        fvm.center[i].dx.push_back(abs(x - x_center));

                        //For y : 
                        int idx_1 = fvm.center[i].shared_face[n][0], idx_2 = fvm.center[i].shared_face[n][1];
                        double y1 = vertex_coordinate_y[idx_1], y2 = vertex_coordinate_y[idx_2];
                        fvm.center[i].dy.push_back(abs(y1 - y2));                  
                    }
                    else {
                        //For x : 
                        int idx_1 = fvm.center[i].shared_face[n][0], idx_2 = fvm.center[i].shared_face[n][1];
                        double x1 = vertex_coordinate_x[idx_1], x2 = vertex_coordinate_x[idx_2];
                        fvm.center[i].dx.push_back(abs(x1 - x2));    
                        
                        //For y : 
                        double y_center = fvm.center[i].center_coordinate[1];
                        double y = fvm.center[neighbor_number].center_coordinate[1];
                        fvm.center[i].dy.push_back(abs(y - y_center));
                        
                    }
                }
                
            }
        }
    }
#endif