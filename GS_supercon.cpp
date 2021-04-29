//======================================================================================================
//
//   Description: A C++ program to solve the Grad-Shafranov equation to find equilibria models
//                with superconducting core
//                                  
//   
//   Version: 1.0
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
//=====================================================================================================



#include <iostream>
#include <math.h>
#include <typeinfo>
#include <fstream>
using namespace std;


int main() {

    int Nr = 51;
    int Nu = 81;
    double u[Nu], r[Nr], theta[Nu];
    double du=2.0/(Nu-1);
    double dr=2.0/(Nr-1);
    int i,j;
    double A[Nr][Nu], B[Nr][Nu], Ac[Nr][Nu], S_nm[Nr][Nu], S_sc[Nr][Nu], Br[Nr][Nu], Bth[Nr][Nu], Bphi[Nr][Nu], Pif[Nr][Nu], Bmag[Nr][Nu], Bcc[Nr][Nu], y[Nr][Nu], Mn[Nr][Nu], Msc[Nr][Nu], fsc[Nr][Nu], dPdA_term[Nr][Nu], fsc_prime[Nr][Nu],yprime[Nr][Nu];
    int error = 1.0;
    double strength;
    double power;
    double den,S2,dS2,Q,w,max1,min1,counter,ncounter,threshold,dth,c0,c1,c2,Dp,dAdr,dAdu;
    int rid,uid,r_cc_id,u_mid_id;
    ofstream outputfile1,outputfile2;
    
    //required input parameters
    threshold = 1e-5;
    power = 1.1;
    strength = 5;
    
    
    //extra parameters
    counter=0;
    ncounter = 1.0;
    S2 = 0;
    dS2 = 0;
      
    Dp = 1.0;
    
    //under-relaxation parameter
    w = 0.1;
    
    // initialize grid
    for (j=0; j<Nu; j++){ 
        u[j] = -1.0 + j*du;
        }
        
    for (i=0; i<Nr; i++){ 
        r[i] = i*dr;
        }
    
        
    //find stellar radius
    for (i=0; i<Nr; i++){
       if (r[i]>1.0){
         rid = i;
         //cout<<"radius is located at="<<i;
         break;
       }  
    }

    //find crust-core interface
    for (i=0; i<Nr; i++){
       if (r[i]>0.9){
         r_cc_id = i;
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

    for (j=0; j<Nu; j++){
       if (u[j]>0.52){
         u_mid_id = j;
         //cout<<"equator is located="<<j;
         break;
       }  
    }
     
    // initialize source
    for (i=0; i<Nr; i++){ 
        for(j=0;j<Nu;j++){
            if (r[i]<1 && r[i]>0.9){
                S_nm[i][j] = pow(r[i],2.0)*(1-pow(u[j],2))*(1-r[i]*r[i]);
                }
            else{
                S_nm[i][j] = 0.0;
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
      
    while (counter<ncounter){     
        
    std::copy(&A[0][0], &A[0][0]+Nr*Nu,&Ac[0][0]);   

     for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (A[i][j]>A[rid][uid]){B[i][j] = strength*pow((A[i][j]-A[rid][uid]),power);}  
            else{ B[i][j] = 0.0 ;}
          
       }
    }

    // calculate magnetic field components
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            theta[j] = acos(u[j]);
            dth = theta[2]-theta[1];
            Br[i][j] = 1.0/r[i]/r[i]/sin(theta[j])*(A[i][j+1]-A[i][j-1])/2/dth;
            Bth[i][j] = -1.0/r[i]/sin(theta[j])*(A[i+1][j]-A[i-1][j])/2/dr;
            Bphi[i][j] = 1.0/r[i]/sin(theta[j])*(B[i][j]);
            Bmag[i][j] = sqrt(Br[i][j]*Br[i][j] + Bth[i][j]*Bth[i][j] + Bphi[i][j]*Bphi[i][j]);
        }
    }
    
    // calculate the different terms for the superconducting core
    c0 = Bmag[r_cc_id][1];
    //cout<<"c0="<<c0<<"\t";
    c1 = (Bmag[r_cc_id][uid]-c0)/A[r_cc_id][uid];
    //cout<<"c1="<<c1<<"\t";
    if ((A[r_cc_id][0]-A[r_cc_id][uid])!=0){
    c2 = (Bmag[r_cc_id][u_mid_id]-c0-c1*A[r_cc_id][u_mid_id])/A[r_cc_id][u_mid_id]/(A[r_cc_id][0]-A[r_cc_id][uid]);}
    else{c2=0;}
    //cout<<"c2="<<c2<<"\t";
    

    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            Bcc[i][j] = c0 + c1*A[i][j] + c2*A[i][j]*(A[i][j]-A[r_cc_id][uid]);
            //cout<<Bcc[i][j]<<" ";
            Mn[i][j] = A[i][j];
            y[i][j] = Dp*Bcc[i][j] + Mn[i][j]; //*4*3.14159;
            //cout<<y[i][j]<<" ";
            Msc[i][j] = (y[i][j] - Bmag[i][j]*Dp);// /4.0/3.14159;
            //cout<<Msc[i][j]<<" ";
            fsc[i][j] = (1.0-r[r_cc_id]*r[r_cc_id])*B[i][j]/Bcc[i][j];
            //cout<<fsc[i][j]<<" ";
            }
    }

    //calculate the Pi factor
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (i<rid-1){
            Pif[i][j] = Bmag[i][j]/(1-r[i]*r[i]);}
            else{ Pif[i][j] = 0.0;}
        }
    }
    


    //calculate the derivative gradPi dot grad alpha
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (i<r_cc_id-1){
            dPdA_term[i][j] = ((Pif[i+1][j]-Pif[i-1][j])/dr*(A[i+1][j]-A[i-1][j])/dr + (Pif[i][j+1]-Pif[i][j-1])/du*(A[i][j+1]-A[i][j-1])/du)/Pif[i][j];}
            else{ dPdA_term[i][j]=0.0;}
            //cout<<dPdA_term[i][j]<<" ";
        }
    }

     //calculate the derivatives of functions present in SC core
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
           dAdr = (A[i+1][j]-A[i-1][j])/dr;
           dAdu = (A[i][j+1]-A[i][j-1])/du;

           if (dAdr==0 || dAdu==0){
                         yprime[i][j] = 0.0;
                        fsc_prime[i][j] = 0.0;
                        }
           else{
           yprime[i][j] = (y[i+1][j]-y[i-1][j])/dr/dAdr + (y[i][j+1]-y[i][j-1])/du/dAdu;          
           fsc_prime[i][j] = (fsc[i+1][j]-fsc[i-1][j])/dr/dAdr + (fsc[i][j+1]-fsc[i][j-1])/du/dAdu;}
          }
       }

   
    //calculate the Source for Superconducting core
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (i<rid-1){
                S_sc[i][j] = (1-r[i]*r[i])*Pif[i][j]*yprime[i][j] + Pif[i][j]*Pif[i][j]*fsc[i][j]*fsc_prime[i][j] - dPdA_term[i][j];     
                //S_sc[i][j] = 0.08*(Pif[i][j]*Pif[i][j]*fsc[i][j]*fsc_prime[i][j] - dPdA_term[i][j]);     
            }
            else{ 
                 S_sc[i][j] = 0.0;}
        }
    }

    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (std::abs(S_sc[i][j])>100.0){S_sc[i][j]=0.0;}
            //cout<<S_sc[i][j]<<" ";
        }
    }
     
    // main loop to update A (i.e. alpha) 
    
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
                
            if (Ac[i][j]>Ac[rid][uid]){
                S2 = Pif[i][j]*Pif[i][j]*strength*strength*power*pow((Ac[i][j]-Ac[rid][uid]),(2.0*power-1.0));  
                dS2 = Pif[i][j]*Pif[i][j]*strength*strength*power*(2*power-1)*pow((Ac[i][j] - Ac[rid][uid]),(2.0*power-2.0));
                
            }
            else{
                S2 = 0.0;
                dS2 = 0.0;
            }
            max1 = std::max(0.0,dS2);
            min1 = std::min(0.0,dS2);
            Q = (S2-dS2*Ac[i][j]) + max1*Ac[i][j];
            den = (2.0/dr/dr + 2.0*(1-u[j]*u[j])/r[i]/r[i]/du/du) - min1;
            A[i][j] = (1-w)*Ac[i][j] + w*((A[i+1][j]+A[i-1][j])/dr/dr + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/r[i]/r[i]/du/du + Q + S_nm[i][j] + 0.01*S_sc[i][j])/den;
                      
       }
    }
  
        
    // calculate the error

    double e = 0;
    for (int i=1; i< Nr-1; i++) {
      for (int j=1; j< Nu-1; j++){
         e += pow((A[i][j]-Ac[i][j]),2)/pow(Ac[i][j],2);
      }
    }

    cout<<"\n error="<<sqrt(e)<<"\n"; 
    
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
     
    outputfile2.open("A_"+std::to_string(int(strength))+ "_supercon.txt");
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
