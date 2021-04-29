//====================================================================================
//
//   Description: A C++ program to solve the Grad-Shafranov equation to find equilibria models
//                                  
//   
//   Version: 4.0
//   Date: Nov 2020
//  
//
//   Equation: GS A = -r**2*sin(theta)*(electron denisty) + Beta*d(Beta)/dA
//               B = strength * (A - A(R,theta))**(power)
//              
//   The code requires you to put an input strength, power, and threshold error
//          
//   Author: Ankan Sur
//   Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
//
//====================================================================================



#include <iostream>
#include <math.h>
#include <typeinfo>
#include <fstream>
using namespace std;


int main() {

    int Nr = 201;
    int Nu = 201;
    double u[Nu], r[Nr];
    double du=2.0/(Nu-1);
    double dr=2.0/(Nr-1);
    int i,j;
    double A[Nr][Nu], Ac[Nr][Nu], S[Nr][Nu];
    int error = 1.0;
    double strength;
    double power;
    double den,S2,dS2,Q,w,max1,min1,counter,threshold;
    //double ,err_arr[1000000],counter_arr[10000000];
    int rid,uid;
    //string error_conv;
    //const char* yes = "yes";
    //const char* no = "no";
    ofstream outputfile1,outputfile2;
    
    //required input parameters
    threshold = 1e-6;
    power = 1.1;
    cout<<"enter the strength=";
    cin>>strength;

    double k;
    cout<<"enter k=";
    cin>>k;

       
    //extra parameters
    counter=0;
    S2 = 0;
    dS2 = 0;
    
    //under-relaxation parameter
    w = 0.9;
    
    // initialize grid
    for (j=0; j<Nu; j++){ 
        u[j] = -1.0 + j*du;
        }
        
    for (i=0; i<Nr; i++){ 
        r[i] = i*dr;
        }
    
    r[0] = 0.001;
        
    //find stellar radius
    for (i=0; i<Nr; i++){
       if (r[i]>1.0){
         rid = i;
         //cout<<"radius is located at="<<i;
         break;
       }  
    }
    //find equator
    for (j=0; j<Nu; j++){
       if (u[j]>0.0){
         uid = j;
         //cout<<"equator is located="<<j;
         break;
       }  
    }
     
    // initialize source
    for (i=0; i<Nr; i++){ 
        for(j=0;j<Nu;j++){
            if (r[i]<1){
                S[i][j] = pow(r[i],2.0)*(1-pow(u[j],2))*(sin(3.14159*r[i])/(3.14159/r[i]))*k;
                }
        }
    }               
    
    // initiliaze A
    for (i = 0; i < Nr; i++) {
        for (j=0; j < Nu; j++){
            A[i][j] = (1-pow(u[i],2));
       }
   }
   
   // set boundaries
    for (i=0;i<Nr;i++){
         A[i][0] = 0.0;
         A[i][Nu-1] = 0.0;
         }
    for (j=0;j<Nu;j++){
        A[0][j] = 0.0;
        A[Nr-1][j] = 0.0;
    }
      
    while (counter<100000000.0){     
        
    std::copy(&A[0][0], &A[0][0]+Nr*Nu,&Ac[0][0]);   
          
    
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
                
            if (Ac[i][j]>Ac[rid][uid]){
                S2 = strength*strength*power*pow((Ac[i][j]-Ac[rid][uid]),(2.0*power-1.0));  
                dS2 = strength*strength*power*(2*power-1)*pow((Ac[i][j] - Ac[rid][uid]),(2.0*power-2.0));
                
            }
            else{
                S2 = 0.0;
                dS2 = 0.0;
            }
            max1 = std::max(0.0,dS2);
            min1 = std::min(0.0,dS2);
            Q = (S2-dS2*Ac[i][j]) + max1*Ac[i][j];
            den = (2.0/dr/dr + 2.0*(1-u[j]*u[j])/r[i]/r[i]/du/du) - min1;
            A[i][j] = (1-w)*Ac[i][j] + w*((A[i+1][j]+A[i-1][j])/dr/dr + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/r[i]/r[i]/du/du + Q + S[i][j])/den;
                      
       }
    }
    
           
    double e = 0;
    for (int i=1; i< Nr-1; i++) {
      for (int j=1; j< Nu-1; j++){
         //cout<<e;
         e += pow((A[i][j]-Ac[i][j]),2)/pow(Ac[i][j],2);
      }
    }

    cout<<sqrt(e)<<"\n"; 
    //err_arr[counter] = sqrt(e);
    //counter_arr[counter] = counter;

    if (sqrt(e)<threshold){
        std::cout<<"Convergence reached: exiting computations"<<"\n";
    	break;
    	}
    	
    if (sqrt(e)>100.0 || isnan(sqrt(e))){
        std::cout<<"Exiting loop, solution diverged"<<"\n";
        break;
        }
          
    
    counter+=1;    
    }          

    
    /*if (error_conv=="yes"){
       outputfile1.open("error.txt");
       for (i=0;i<counter;i++){
        std:cout<<err_arr[i]<<"\t"<<counter_arr[i]<<"\n";    
        }
       outputfile1.close();
     }*/
    
    
    outputfile2.open("A_"+std::to_string(int(strength))+ "_rho_n1.txt");
    cout<<"Produced output file"<<"\n";

    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile2<<A[count][index]<<" ";  
            outputfile2<<endl;                          
        }
    outputfile2.close(); 
    
    
    return 0;
}
