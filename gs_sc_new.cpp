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
    int Nu = 51;
    double u[Nu], r[Nr], theta[Nu];
    double du=2.0/(Nu-1);
    double dr=2.0/(Nr-1);
    int i,j;
    double A[Nr][Nu], B[Nr][Nu], Ac[Nr][Nu], S[Nr][Nu], S_nm[Nr][Nu], S_sc[Nr][Nu], Br[Nr][Nu], Bth[Nr][Nu], Bphi[Nr][Nu], Pif[Nr][Nu], Bmag[Nr][Nu], Bcc[Nr][Nu], y[Nr][Nu], Mn[Nr][Nu], Msc[Nr][Nu], fsc[Nr][Nu], dPdA_term[Nr][Nu], fsc_prime[Nr][Nu],yprime[Nr][Nu];
    double err_arr[Nr][Nu], res[Nr][Nu], term1[Nr][Nu], term2[Nr][Nu], term3[Nr][Nu];
    double e = 1.0;
    double strength;
    double power;
    double den,S2,dS2,Q,w,max1,min1,threshold,dth,c0,c1,c2,Dp,dAdr,dAdu;
    int rid,uid,r_cc_id,u_mid_id;
    ofstream outputfile1,outputfile2;
    
    //required input parameters
    threshold = 1e-5;
    power = 1.1;
    strength = 5;
    
    double h_c=0.9;
    //extra parameters
    int counter=0;
    int ncounter = 10;
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
         rid = i-1;
         break;
       }  
    }
    cout<<"radius ar r="<<r[rid]<<"\n";

    //find crust-core interface
    for (i=0; i<Nr; i++){
       if (r[i]>0.9){
         r_cc_id = i-1;
         break;
       }  
    }


    cout<<"Crust core boundary at r="<<r[r_cc_id]<<"\n";
     
    //find equator
    for (j=0; j<Nu; j++){
       if (u[j]>0.0){
         uid = j-1;
         //cout<<"equator is located="<<j;
         break;
       }  
    }

   for (j=0; j<Nu; j++){
       if (u[j]>0.78){
         u_mid_id = j-1;
         break;
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
        A[Nr-1][j] = (1-u[j]*u[j])/r[Nr-1];
    }
      
    while (e>threshold){     
    //while (counter<ncounter){   
    std::copy(&A[0][0], &A[0][0]+Nr*Nu,&Ac[0][0]);   


    for (i=0; i<Nr; i++){ 
        for(j=0;j<Nu;j++){
            if (Ac[i][j]>Ac[rid][uid]){
                B[i][j] = strength*pow((A[i][j]-A[rid][uid]),power);
                }
            else{
                B[i][j] = 0.0 ;
                } 
        }
    }


     // calculate magnetic field components
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            dr = r[i]-r[i-1];
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
    c1 = (Bmag[r_cc_id][uid]-c0)/A[r_cc_id][uid];
    if ((A[r_cc_id][0]-A[r_cc_id][uid])!=0){
    c2 = (Bmag[r_cc_id][u_mid_id]-c0-c1*A[r_cc_id][u_mid_id])/A[r_cc_id][u_mid_id]/(A[r_cc_id][u_mid_id]-A[r_cc_id][uid]);}
    else{c2=0;}
    
     //calculate the Pi factor
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            if (i<=r_cc_id){
            Pif[i][j] = Bmag[i][j]/(1-r[i]*r[i])/h_c;}
            else{ Pif[i][j] = 0.0;}
        }
    }

     for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
          if (i<=r_cc_id){
         Bcc[i][j] = c0 + c1*A[i][j] + c2*A[i][j]*(A[i][j]-A[r_cc_id][uid]);
         }
       else{
       Bcc[i][j] = 0.0;}
   }
  }


    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            Mn[i][j] = A[i][j];
            //y[i][j] = Bcc[i][j] + Mn[i][j]/h_c; 
            if (i<=r_cc_id){
               y[i][j] = Bcc[i][j] + Mn[i][j]/h_c; 
               Msc[i][j] = (y[i][j] - Bmag[i][j]*Dp);            
               if (Pif[i][j]!=0){
                  dPdA_term[i][j] = ((Pif[i+1][j]-Pif[i-1][j])*(A[i+1][j]-A[i-1][j])/dr/dr + (1-u[j]*u[j])*(Pif[i][j+1]-Pif[i][j-1])*(A[i][j+1]-A[i][j-1])/du/du/r[i]/r[i])/Pif[i][j]/4.0;
                  //fsc[i][j] = r[i]*sqrt(1-u[j]*u[j])*B[i][j]/Pif[i][j];
                  //fsc[i][j] = h_c*(1-r[r_cc_id]*r[r_cc_id])*B[i][j]/Bcc[i][j];
               }
               else if (Bcc[i][j]!=0){fsc[i][j] = h_c*(1-r[r_cc_id]*r[r_cc_id])*B[i][j]/Bcc[i][j];}
               
               else{
                   fsc[i][j]=0.0;
                   dPdA_term[i][j]=0.0;
                   }
            }
            else{
               y[i][j] = 0.0;
               fsc[i][j]=0.0;
               Msc[i][j] = 0.0;
               dPdA_term[i][j]=0.0;
            }
        }
   }

  
   for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){


           dAdr = (Ac[i+1][j]-Ac[i-1][j])/dr/2.0;
           dAdu = (Ac[i][j+1]-Ac[i][j-1])/du/2.0;         
    
           

            if (i<=rid && i>r_cc_id){


                if (dAdr==0 || dAdu==0){
                S_nm[i][j] = r[i]*r[i]*(1-u[j]*u[j])*(1-r[i]*r[i]);}
                
                else{
                S_nm[i][j] = r[i]*r[i]*(1-u[j]*u[j])*(1-r[i]*r[i]) + B[i][j]*((B[i+1][j]-B[i-1][j])/dr/dAdr + (B[i][j+1]-B[i][j-1])/du/dAdu);
                }
            }
            else{
                S_nm[i][j] = 0.0;
            }


            if (i<=r_cc_id){


                 if (dAdr==0 || dAdu==0){

                        yprime[i][j] = 0.0;
                        fsc_prime[i][j] = 0.0;
                        }
                 else{
                        yprime[i][j] = c1 + c2*(2*Ac[i][j]-Ac[r_cc_id][u_mid_id]) + 1.0/h_c;
                        fsc_prime[i][j] = (fsc[i+1][j]-fsc[i-1][j])/(2.0*dr)/dAdr + (fsc[i][j+1]-fsc[i][j-1])/(2.0*du)/dAdu;
                  }




                term1[i][j] = r[i]*r[i]*(1-u[j]*u[j])*(1-r[i]*r[i])*Pif[i][j]*yprime[i][j];
                term2[i][j] = Pif[i][j]*Pif[i][j]*fsc[i][j]*fsc_prime[i][j];
                term3[i][j] = dPdA_term[i][j];
                S_sc[i][j] = term1[i][j]+term2[i][j]-term3[i][j];     
               
            }
            else{ 
                 term1[i][j] = 0.0;
                 term2[i][j] = 0.0;
                 term3[i][j] = 0.0;
                 S_sc[i][j] = 0.0;}
            
            if (S_sc[i][j] != S_sc[i][j]) {S_sc[i][j] = 0.0;}
            if (std::isinf(S_sc[i][j])){S_sc[i][j] = 0.0;}


            S[i][j] = S_sc[i][j] + S_nm[i][j];  

    }
   }

   /*for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            S[i][j] = (S_sc[i+1][j] + S_sc[i-1][j] + S_sc[i-1][j-1]  + S_sc[i-1][j+1] + S_sc[i][j-1] + S_sc[i][j+1] + S_sc[i+1][j-1] + S_sc[i+1][j+1]+S_sc[i][j])/9.0;
    }
   }*/

      /*for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){


           dAdr = (Ac[i+1][j]-Ac[i-1][j])/dr/2.0;
           dAdu = (Ac[i][j+1]-Ac[i][j-1])/du/2.0;         
    
                     
           if (r[i]<1 && r[i]>0.0){
                if (dAdr==0 || dAdu==0){
                S_nm[i][j] = r[i]*r[i]*(1-u[j]*u[j])*(1-r[i]*r[i]);}

            else{
                S_nm[i][j] = r[i]*r[i]*(1-u[j]*u[j])*(1-r[i]*r[i]) + B[i][j]*((B[i+1][j]-B[i-1][j])/dr/dAdr + (B[i][j+1]-B[i][j-1])/du/dAdu);
                }
            }
            else{
                S_nm[i][j] = 0.0;
            }

    
               }
         }*/
           


    
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){

            den = (2.0/dr/dr + 2.0*(1-u[j]*u[j])/r[i]/r[i]/du/du);
            res[i][j] = ((A[i+1][j]+A[i-1][j])/dr/dr + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/r[i]/r[i]/du/du + S[i][j]);
            A[i][j] = (1-w)*Ac[i][j] + w*((A[i+1][j]+A[i-1][j])/dr/dr + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/r[i]/r[i]/du/du + S[i][j])/den;
                     
       }
    }

    // set boundaries
    for (i=0;i<Nr;i++){
         A[i][0] = 0.0;
         A[i][Nu-1] = 0.0;
         }
    for (j=0;j<Nu;j++){
        A[0][j] = 0.0;
        A[Nr-1][j] = A[rid][uid]*(1-u[j]*u[j])/r[Nr-1];
    }
    
        
    // calculate the error

    double e = 0;
    for (int i=1; i< Nr-1; i++) {
      for (int j=1; j< Nu-1; j++){
         err_arr[i][j] = sqrt(pow((A[i][j]-Ac[i][j]),2)/pow(Ac[i][j],2));
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
          
    /*if (counter%100==0){
    outputfile1.open("res_"+std::to_string(int(counter))+".txt");

    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile1<<err_arr[count][index]<<" ";  
            outputfile1<<endl;                          
        }
    outputfile1.close(); 
    }*/
    

    counter+=1;    
    }      

    outputfile1.open("Pi.txt");

    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile1<<Pif[count][index]<<" ";  
            outputfile1<<endl;                          
        }
    outputfile1.close(); 
           
    
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
