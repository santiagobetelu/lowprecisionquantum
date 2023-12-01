/////////////////////////////////////////////////////////////////////////////
//  TEST PROGRAM FOR LOWPRECISION: RANDOM ALGORITH IN TABLE VII OF PAPER
//  (c) Santiago Ignacio Betelu, Denton 2020
//  mpicc -march=native -Ofast test-smallcomplex.c -o test -lm -Wall
//  mpirun ./test
////////////////////////////////////////////////////////////////////////////////
#define QUBITS 20 

// TRIPLET E,F,A FOR Eq. (3)
#define Esize  4 
#define Fsize  9
#define Asize  11

#include "smallcomplex-library.c"

////////////////////////////////////////

#define CYCLES 4

void qcprogram(){
   int8_t sign[3*QUBITS*CYCLES]; // signs for rotations
   int64_t i,z,x,x0,q1;
   double err,err0,theo,epsilonc2,G,p,mu, entropy0,entropy,theta,lambda,phi;
   complex double c0,c1;

   // create a random list of signs for the rotation U3 (see Table VII)
   srand(time(NULL)/60); // VERY IMPORTANT ALL NODES HAVE SAME SEED
   srand48(time(NULL)/60+inode); // DIFFERENT SEED PER NODE, THIS IS FOR NORMALIZATION
   for(i=0;i<3*CYCLES*QUBITS;i++) sign[i]= 1-2*(rand()%2);
   
   // initial condition Eq. (before 28)
   x0= rand()%N;
   for(z=0;z<N/nnodes;z++){
       x= z+inode*(N/nnodes);
       c0=0.0;
       if( x==x0 ) c0=1.0;
       amplitudeset(c ,z, c0 );
   }

   G=0; // number of noisy gates
   for(i=0; i<CYCLES; i++){
        for(q1=0;q1<QUBITS;q1++){

           theta=  pi/2*sign[3*(i*QUBITS+q1)+0];  
           lambda= pi/4*sign[3*(i*QUBITS+q1)+1]; 
           phi=    pi/4*sign[3*(i*QUBITS+q1)+2]; 
 
           U3(q1, theta ,lambda, phi );     
           CNOT(q1, (q1+1)%QUBITS );
           G=G+1;     
           if(inode==0) printf("%ld ",q1);  
        }
        if(inode==0) printf(" cycle:%ld\n",i);
   }          
   normalize();
 
   // calculate entropy
   entropy0=0.0; 
   for(z=0;z<N/nnodes;z++){
      c0= amplitude(c,z);
      p= creal( c0*conj(c0) );
      if(p>0.0) entropy0= entropy0- p*log(p);
   }
   MPI_Allreduce(&entropy0,&entropy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if(inode==0) printf("S:%g exact:%g\n",entropy, log(1.0*N)-1.0+0.5772156649 );

   for(i=CYCLES-1; i>=0; i--){
        for(q1=QUBITS-1;q1>=0;q1--){

           theta= pi/2*sign[3*(i*QUBITS+q1)+0];  
           lambda=pi/4*sign[3*(i*QUBITS+q1)+1];
           phi=   pi/4*sign[3*(i*QUBITS+q1)+2];
 
           CNOT( q1, (q1+1)%QUBITS );
           U3(q1, -theta ,-phi, -lambda );     
           G=G+1;    
           if(inode==0) printf("%ld ",q1);
        }
        if(inode==0) printf(" cycle:%ld\n",i);
   }
   normalize();

   // verify result with Eqs. (28,15)
   err0=0.0; 
   for(z=0;z<N/nnodes;z++){
      c1=0.0;
      if( z+inode*(N/nnodes)==x0 ) c1=1.0;
      c0= amplitude(c,z);
      err0= err0+ creal((c0-c1)*conj(c0-c1));
   }
   MPI_Allreduce(&err0,&err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   if(inode==0){
            mu= exp(-(pow(2.0, Esize)-pow(0.5, Fsize)));  // smallest representable complex number Eq. (5)
            phi= 1.0-(mu*mu*N+1.0)*exp(-mu*mu*N); // loss of normalization error
            epsilonc2=  phi+ (1.0-phi)*(pow(0.5,2.0*Fsize)+4*pi*pi*pow(0.5,2.0*Asize))/12.0; // Eq. (10)
            theo= G*epsilonc2; // Eq. (15)
            printf("G:%g   sigma^2_theo:%g   sigma^2_actual:%g\n", G, theo , err );
            printf("If all is OK the errors should be small, but not necesarily equal because this run is non-random\n");
   }
   return;
}
////////////////////////////////////////////////////////////////////////////////
