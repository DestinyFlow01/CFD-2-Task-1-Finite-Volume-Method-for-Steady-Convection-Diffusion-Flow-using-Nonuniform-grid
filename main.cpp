#include<iostream>
#include<vector>
#include<string>
#include "FVM.h"
using namespace std;

int main() {
    PrintInfo();

    //Defining variable for simulation :
    FVM FV_case;
    
    //Initialization
    FV_case.Initialization();

    #ifdef General_Method
        Find_Neighbor(FV_case);
    #endif

    cout<<"Database for the generated mesh : \n";
    PrintCenter(FV_case);
    cout<<"---------------------------------------------------------------------------------------------------------------\n";

    //Finite volume method : 
    FV_case.Finite_Volume_Method();

    //Output Result : 
    OutputCSV(FV_case, "Output Result.csv");

}