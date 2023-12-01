/////////////////////////////////////////////////////////////////////////////
//  LIBRARY FILE
//  LOW PRECISION LIBRARY FOR QUANTUM CIRCUIT SIMULATION
//  8-bit aligned version
//  (c) Santiago Ignacio Betelu, Denton 2020
//    __   ____  _      __ ___   ___   ____ _____ ____ ____ ____ ____   _  __
//   / /  / __ \| | /| / // _ \ / _ \ / __// ___//  _// __//  _// __ \ / |/ /
//  / /__/ /_/ /| |/ |/ // ___// , _// _/ / /__ _/ / _\ \ _/ / / /_/ //    / 
// /____/\____/ |__/|__//_/   /_/|_|/___/ \___//___//___//___/ \____//_/|_/  
//                        ____   __  __ ___    _  __ ______ __  __ __  ___  
//                       / __ \ / / / // _ |  / |/ //_  __// / / //  |/  /  
//                      / /_/ // /_/ // __ | /    /  / /  / /_/ // /|_/ /   
//                      \___\_\\____//_/ |_|/_/|_/  /_/   \____//_/  /_/    
//  
// Many thanks to Datavortex Technologies that supported this work and 
// provided the computer system Hypatia used in most of the development and 
// to the Texas Advanced Computing Center (TACC) for providing 
// access to Stampede2.                                                                                        
//////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <mpi.h>
#include <time.h>

#ifndef QUBITS 
   // the default if not defined in driver program "qcprogram"
   #define QUBITS 20ul 
#endif

#define N (1ul<<QUBITS)

#define pi 3.14159265358979323846l

void qcprogram();

uint8_t *c, *buffer; // quantum amplitudes
int64_t BUFFERSIZE,NBUFFERS,NODEBITS,nnodes,inode;
double globalnorminv=1.0; // for renormalization

//////////////////////////////////////////////
// z= x+I*y,  log z= log r+ I*arg(z) 
// where log r<0 thus store as 
//     log r= -(R+F/2^LOGFRAC)  and  arg(z)= pi + pi T/2^Asize
//     d= [R.F][A] = [Esize.Fsize] [Asize]

// TRIPLET E,F,A FOR Eq. (3)
#ifndef Esize
   #define Esize  4 
   #define Fsize  9
   #define Asize  11
#endif

#define COMPLEXBITS (Esize+Fsize+Asize)                 // bits per complex
#define COMPLEXBYTES ( (COMPLEXBITS/8) + ((COMPLEXBITS%8)>0) ) // bytes per complex
#define ZEROCODE (0xfffffffffffffffful>>(64-Esize-Fsize-Asize))
#define MASKARG ((1ul<<Asize)-1ul)

/////////////////////////////////////////////////
// UNOPTIMIZED VERSION OF THE CONVERSION FROM DOUBLE COMPLEX TO LOWPRECISION (EQ. (4) IN PAPER)
uint64_t complex2compact0(complex double z){ 
   double x,y,r2,logr,q;
   uint64_t d,logri, argi;

   x=creal(z); y=cimag(z);
   r2= x*x+y*y;
   if(r2==0.0) return(ZEROCODE);
   logr=  0.5*log(r2) ;
   logri= round(-logr*( 1<<Fsize ) ); // rounded and non-negative
   if(r2>1.0) logri=0;  // correct errors that make amplitude larger than one
   if( logri > ((1<<(Esize+Fsize))-1) ) return(ZEROCODE);
   q= carg(z);
   if(q<0.0 ) q=q+2.0*pi;
   if(q>2*pi) q=q-2.0*pi;
   argi= round( q/(2.0*pi)*( 1<<Asize ));    // rounded and non-negative
   argi= argi & MASKARG;  // periodicity
   d= logri;
   d= (d<<Asize) | argi;
   return(d);
}
///////////////////////////////////////////////
// UNOPTIMIZED VERSION OF THE CONVERSION FROM LOWPRECISION TO DOUBLE COMPLEX Eq. (3)
complex double compact2complex0(uint64_t d){ // unoptimized
   double r,q;
   complex double z;
   uint64_t argi, logri;
   if( d==ZEROCODE ) return(0.0);
   logri= d>>Asize;
   argi = d & MASKARG;
   r= exp( -1.0*logri/(1<<Fsize) ); 
   q= 2.0*pi*argi/( 1<<Asize ) ;
   z= r*cexp(I*q);  
   return(z);
}
/////////////////////////////////////////////////////////////////
// MY FIRST ATTEMPT AT OPTIMIZATION, VERY CRUDE AND NOT SO FAST...
uint64_t complex2compact(complex double z){
   #define NEsizeERP 1024 
   #define NATANINTERP 1024
   #define NATAN (1<<Asize)            // must be a power of 2
   #define NLOGR (1<<(Fsize+Esize)) // must be a power of 2
   float x,y,v,dv,logr,r2; // float for speed
   static float logrt[NEsizeERP+2], atant[NATANINTERP+2];
   static int first=1;
   uint64_t d;
   int32_t i,j,s,logri,argi,m,e; // signed

   if( first==1 ){
       first=0; // fill lookup tables at the beginning
       for(i=0; i<=NEsizeERP+1; i++){
          logrt[i]= log(1.0+1.0*i/NEsizeERP); // log(1+mantissa) [0,1]
       }
       for(i=0; i<=NATANINTERP+1; i++){
          atant[i]= atan( 1.0*i/NATANINTERP )*4.0/pi *(NATAN/8); // arctan interpolation table argument 0<=t<=1, range 0<=q<=pi/4
       }
   }
   x=crealf(z); y=cimagf(z);
   r2= x*x+y*y;
   if(r2==0.0) return(ZEROCODE);
   if( r2>=1.0 ) {
      logri=0;
   }else{
      memcpy(&s,&r2,4);
      e= (s>>23)-127; // exponent, the sign always 0
      m= s & 0x7FFFFF; // mantissa
      v=  1.0*NEsizeERP*m/0x7FFFFF; // fractionary mantissa 0<=v<1
      j=  v;    // table index for the mantissa
      dv= v-j;  // fractionary difference
      logr= (1.0f-dv)*logrt[j] + dv*logrt[j+1] +log(2.0)*e  ; 
      logri= round(  -0.5f*(1<<Fsize)*logr  )  ; 
      if( logri >= NLOGR ) return(ZEROCODE);
   }
   if( fabsf(x)<fabsf(y) ){
       v=  fabsf(x/y);
       j=  v*NATANINTERP;
       dv= v*NATANINTERP-j;
       argi= roundf( (1.0f-dv)*atant[j] + dv*atant[j+1] );
       argi= NATAN/4-argi;
   }else{
       v=  fabsf(y/x);
       j=  v*NATANINTERP;
       dv= v*NATANINTERP-j;
       argi= roundf( (1.0f-dv)*atant[j] + dv*atant[j+1] );
   }
   if( x<0.0f )
      if( y<0.0f )
         argi= NATAN/2+argi;
      else
         argi= NATAN/2-argi;
   else
      if( y<0.0f ) argi= NATAN-argi;
   argi= argi & MASKARG;  // periodicity
   d= logri;
   d= (d<<Asize) | argi;
   return(d);
}
/////////////////////////////////////////////////////////////////
// OPTIMIZED VERSION FROM LOWPRECICION TO COMPLEX, WORKS GREAT
complex double compact2complex(uint64_t d){
   static int first=1;
   static complex float argtable[NATAN];
   static float rhotable[NLOGR];
   double r,q; // for precomputing only
   complex double z;
   uint32_t argi, logri;
   if(first==1){ // precompute rho and theta for speed
      first=0;
      for(argi=0; argi<NATAN; argi++){
            q= 2.0*pi*argi/( 1<<Asize );
            argtable[argi]= cexp(I*q);
      }
      for(logri=0; logri<NLOGR; logri++){
            r= exp( -1.0*logri/(1<<Fsize) );
            rhotable[logri]= r;
      }
   }
   d= d&ZEROCODE;
   if( d==ZEROCODE ) return(0.0);
   logri= d>>Asize;
   argi = d & MASKARG;
   z= rhotable[logri]*argtable[argi];  // can be done with lookup table
   return(z);
}
///////////////////////////////////////////////////////
// RETURNS THE COMPLEX DOUBLE STORED IN ca[i] WITH LOWPRECISION
complex double amplitude(uint8_t *ca, int64_t i){ // obtains complex amplitude c0=c[i]
  complex double c0=0.0;
  uint64_t d,k;
  if( i<N/nnodes ){
     k= i*COMPLEXBYTES;
     d=0;
     memcpy( &d, ca+k, COMPLEXBYTES);
     c0=compact2complex( d );
  }
  return(c0);
}
////////////////////////////////////////////////////////
void amplitudeset( uint8_t *ca, int64_t i, complex double c0 ){ // sets complex amplitude c[i]=c0
  uint64_t d,k;
  if( i<N/nnodes ){
     k= i*COMPLEXBYTES;
     d= complex2compact(c0);
     memcpy( ca+k, &d, COMPLEXBYTES);
  }
  return;
}
///////////////////////////////////////////////////////
// EFFICIENT COPY c1[i] TO c2[j] IN LOWPRECISION
void amplitudecopy( uint8_t *c1, int64_t i, uint8_t *c2, int64_t j ){ 
  if( i<N/nnodes && j<N/nnodes )
      memcpy( c1+(i*COMPLEXBYTES), c2+(j*COMPLEXBYTES), COMPLEXBYTES);
  return;
}
//////////////////////////////////////////////////////////
// NORMALIZATION PROCEDURE FOR AMPLITUDES SECTION IV-B IN PAPER Eq. (23)
void normalize(){
   int64_t z;
   double p,p0,delta;
   complex double c0;

   p0=0.0;
   for(z=0;z<N/nnodes;z++){
       c0= amplitude(c,z);
       p0= p0+creal(c0*conj(c0));
   }
   MPI_Allreduce(&p0,&p, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   globalnorminv= 1.0/sqrt(p);
 
   for(z=0;z<N/nnodes;z++){
      delta= pow(0.5,1.0*Fsize)*(0.5-drand48());
      c0= amplitude(c,z);
      c0= c0*globalnorminv*( 1.0 +delta*(1.0+0.5*delta*(1.0+delta/6.0) )); // ~ exp(delta)
      amplitudeset(c,z,c0);
   }
   return; 
}
/////////////////////////////////////////////////////////////////////////////
//  Quantum numstates addressing example with 4 nodes.
//  1- The QUBITS-NODEBITS least significant bits can be swapped within each node.
//  2- The NODEBITS most significant digits is node number
//     NODE
//       |  Local bits
//  c0   00 000           c16  10 000
//  c1   00 001           c17  10 001
//  c2   00 010           c18  10 010
//  c3   00 011  N0       c19  10 011  N2
//  c4   00 100           c20  10 100
//  c5   00 101           c21  10 101
//  c6   00 110           c22  10 110
//  c7   00 111           c23  10 111
//       ......                ......
//  c8   01 000           c24  11 000
//  c9   01 001           c25  11 001
//  c10  01 010           c26  11 010
//  c11  01 011  N1       c27  11 011  N3
//  c12  01 100           c28  11 100
//  c13  01 101           c29  11 101
//  c14  01 110           c30  11 110
//  c15  01 111           c31  11 111
//       ......                ......
//////////////////////////////////////////////////////////////////////////////
//  H= | 1  1 |
//     | 1 -1 | /sqrt(2)
void H(uint64_t qubit1){  // Hadamard gate acting on qubit1
    uint64_t x,y,q,chunk,node;
    complex double c0,c1,cx,cy;
    MPI_Status status;
    // compute globalnorminv in the fly
    if(qubit1< QUBITS-NODEBITS){
       for(x=0; x<N/nnodes; x++){
           y= x^(1ul<<qubit1);   // XOR for exchange
           if(x<y){
              cx=  amplitude(c,x);
              cy=  amplitude(c,y); 
              c1=  (cx-cy)*sqrt(0.5);
              c0=  (cx+cy)*sqrt(0.5);
              amplitudeset(c,x,c0);
              amplitudeset(c,y,c1);
           }
       }
    }else{
       node= inode^(1ul<<(qubit1-(QUBITS-NODEBITS))); // which node to communicate
       for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
            MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                            buffer,              BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                          MPI_COMM_WORLD, &status);
            if( inode&(1ul<<(qubit1-(QUBITS-NODEBITS))) ){
                for(q=0; q<BUFFERSIZE; q++){
                   c1= -(amplitude(c,chunk+q)-amplitude(buffer,q))*sqrt(0.5);
                   amplitudeset(c,chunk+q,c1);
                }
            }else{
                for(q=0; q<BUFFERSIZE; q++){
                   c1= (amplitude(c,chunk+q)+amplitude(buffer,q))*sqrt(0.5);
                   amplitudeset(c,chunk+q,c1);
                }
            }
       }
    }
    return; 
}
///////////////////////////////////////////////////////////////////////////////
void U3(int64_t qubit, double theta, double lambda, double phi){  // rotation about angles z in Bloch sphere Eq. (25)
    uint64_t x,y,q,chunk,node;
    double complex explambda,expphi,cx,cy,c0,c1;
    double costheta,sintheta;
    MPI_Status status;
    //
    sintheta=sin(0.5*theta); costheta=cos(0.5*theta);
    explambda= cexp(I*lambda); expphi= cexp(I*phi);
    if(qubit< QUBITS-NODEBITS){
       for(x=0;x<N/nnodes;x++){
           y= x^(1ul<<qubit);   // XOR for the exchange
           if(x<y){
              cx=  amplitude(c,x);
              cy=  amplitude(c,y);
              c0=  cx*costheta-cy*explambda*sintheta;
              c1=  cx*expphi*sintheta+cy*explambda*expphi*costheta;
              amplitudeset(c,x,c0);
              amplitudeset(c,y,c1);
           }
       }
    }else{
       node= inode^(1ul<<(qubit-(QUBITS-NODEBITS))); // which node to communicate
       for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
            MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                            buffer,              BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                          MPI_COMM_WORLD, &status);
            if( inode&(1ul<<(qubit-(QUBITS-NODEBITS))) ){
                for(q=0; q<BUFFERSIZE; q++){
                   c1=  amplitude(buffer,q)*expphi*sintheta +amplitude(c,chunk+q)*explambda*expphi*costheta;
                   amplitudeset(c,chunk+q,c1);
                }
            }else{
                for(q=0; q<BUFFERSIZE; q++){
                   c1= amplitude(c,chunk+q)*costheta-amplitude(buffer,q)*explambda*sintheta;
                   amplitudeset(c,chunk+q,c1);
                }
            }
       }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////
void SWAP(int64_t qubit1, int64_t qubit2){  // SWAP between qubit1 and qubit2, qubit1!=quibit2
    int64_t x,y,b1,b2,chunk,q;
    int node;
    uint8_t aux[COMPLEXBYTES];
    MPI_Status status;
    //
    memset(aux,0,COMPLEXBYTES);
    if(qubit1>qubit2){ // sort qubit1 < qubit2
        q=qubit1;
        qubit1=qubit2;
        qubit2=q;
    }
    if(qubit2<QUBITS-NODEBITS && qubit1<QUBITS-NODEBITS){
        for(x=0; x<N/nnodes; x++){
           b1= (x>>qubit1)&1ll;
           b2= (x>>qubit2)&1ll;
           if(b1!=b2){
              y= (x^(1ll<<qubit1))^(1ll<<qubit2);
              if(y>x){ // to avoid overwriting previously computed
                 amplitudecopy( aux,0, c,x);
                 amplitudecopy( c,x, c, y);
                 amplitudecopy( c,y, aux,0);
              }
           }
        }
    }else if(qubit1 >= QUBITS-NODEBITS && qubit2 >= QUBITS-NODEBITS) { // in this case swap all array alements with another node
        x=  inode*(N/nnodes);
        b1= (x>>qubit1)&1ul;
        b2= (x>>qubit2)&1ul;
        if( b1!=b2 ){
           node= inode^(1ul<<(qubit2-(QUBITS-NODEBITS)));  // here qubit2 >= QBITS-NODEBITS for sure
           node=  node^(1ul<<(qubit1-(QUBITS-NODEBITS)));
           for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
               MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                            buffer,                 BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                               MPI_COMM_WORLD, &status);
               for(q=0; q<BUFFERSIZE; q++){
                   amplitudecopy(c,chunk+q, buffer,q);
               }
           }
        }
   }else{  // qbit1 inside same node but qbit2 in another node
           //fprintf(stderr,"LH ");
           node= inode^(1ul<<(qubit2-(QUBITS-NODEBITS)));  // here qubit2 >= QBITS-NODEBITS for sure
           if(node>=nnodes) fprintf(stderr,"badnode2");
           x= node*(N/nnodes);
           b2= (x>>qubit2)&1ul;
           for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
                   MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                                buffer,                 BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                                  MPI_COMM_WORLD, &status);
                   for(q=0; q<BUFFERSIZE; q=q+1){
                       x= chunk+q; // received register
                       b1= (x>>qubit1)&1ul;
                       y= (chunk+q)^(1ul<<qubit1);  // guaranteed y<x
                       if( b1!=b2 ) amplitudecopy(c,y, buffer,q);
                   }
           }
    }
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void CP(int64_t qubit1, int64_t qubit2, int64_t phaseexp, int64_t sign){  // PHASE between control qubit1 and qubit2, qubit1!=quibit2, phase= pi/2^k
    int64_t x,q,b1,b2; 
    complex double expphase,c1,delta;
    double phase;
    phase= sign*pi*pow(2.0,-1.0*phaseexp);
    expphase= cexp(I*phase);
    for(q=0;q<N/nnodes;q++){
       x= q+inode*(N/nnodes);
       b1= ((x>>qubit1)&1ul);
       b2= ((x>>qubit2)&1ul);
       if( b1 && b2 ){
          c1=amplitude(c,q);
          c1= c1*expphase;
          if( phaseexp>=Asize ){ // correct biased errors
              delta= 2.0*pi*I*pow(0.5,1.0*Asize)*(0.5-drand48()); // small unbiased complex number
              if( phaseexp==Asize ) delta=0.5*delta; // when in the middle of interval
              c1= c1*( 1.0 +delta*(1.0+0.5*delta*(1.0+delta/6) )); // ~ exp(delta)
          }
          amplitudeset(c,q,c1);
       }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void CNOT(uint64_t qubit1, uint64_t qubit2){  // CNOT between control qubit1 and qubit2, qubit1!=quibit2
    uint64_t x,y,b1,q,chunk,node;
    uint8_t aux[COMPLEXBYTES];
    MPI_Status status;
    //
    memset(aux,0,COMPLEXBYTES);
    if(qubit2< QUBITS-NODEBITS){ // acts inside this node
       for(x=0; x<N/nnodes; x++){
           q= x+inode*(N/nnodes);
           y= x^(1ull<<qubit2);   // NOT on bit qubit2
           b1= ((q>>qubit1)&1ul);
           if( b1 && (x<y) ){
              amplitudecopy( aux,0, c,x);
              amplitudecopy( c,x, c, y);
              amplitudecopy( c,y, aux,0);
           }
       }
    }else{
        node= inode^(1ul<<(qubit2-(QUBITS-NODEBITS))); // which node to communicate
        for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
            MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                            buffer,              BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                          MPI_COMM_WORLD, &status);
            for(q=0; q<BUFFERSIZE; q++){
                x= (chunk+q) +inode*(N/nnodes);
                if( (x>>qubit1)&1ul ){
                    amplitudecopy( c, chunk+q,  buffer, q);
                }
            }
        }
    }
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void NOT(uint64_t qubit1){  // CNOT between control qubit1 and qubit2, qubit1!=quibit2
    uint64_t x,y,q,chunk,node;
    uint8_t aux[COMPLEXBYTES];
    MPI_Status status;
    //
    memset(aux,0,COMPLEXBYTES);
    if(qubit1< QUBITS-NODEBITS){ // acts inside this node
       for(x=0; x<N/nnodes; x++){
           y= x^(1ull<<qubit1);   // NOT on bit qubit2
           if( (x<y) ){
              amplitudecopy( aux,0, c,x);
              amplitudecopy( c,x, c, y);
              amplitudecopy( c,y, aux,0);
           }
       }
    }else{
        node= inode^(1ul<<(qubit1-(QUBITS-NODEBITS))); // which node to communicate
        for(chunk=0; chunk<N/nnodes; chunk=chunk+BUFFERSIZE){
            MPI_Sendrecv(&c[chunk*COMPLEXBYTES], BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, node,
                            buffer,              BUFFERSIZE*COMPLEXBYTES, MPI_UINT8_T, node, inode,
                                                                          MPI_COMM_WORLD, &status);
            for(q=0; q<BUFFERSIZE; q++){
                x= (chunk+q) +inode*(N/nnodes);
                amplitudecopy( c, chunk+q,  buffer, q);
            }
        }
    }
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){
   int64_t aux;
   struct timespec tim0,tim1;
   double time;

   if( COMPLEXBITS>50) { printf("COMPLEXBITS>50\n"); exit(1); } // it would not make sense
   if ( QUBITS>63 ) { printf("QUBITS>63\n"); exit(1); }
   clock_gettime(CLOCK_REALTIME,&tim0);
   
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,(int*)&nnodes);
   MPI_Comm_rank(MPI_COMM_WORLD,(int*)&inode);

   NODEBITS=0;
   aux=1;
   while(aux<nnodes){ aux= (aux<<1); NODEBITS= NODEBITS+1; }
   if(aux!=nnodes){
      fprintf(stderr,"ERROR: Number of nodes has to be a power of 2\n");
      exit(1);
   }
   if(inode==0){
       printf("Ranks:%lu  qubits:%d  Bits per complex:%d  E:%d F:%d A:%d  Nodebits:%lu\n", nnodes, QUBITS, COMPLEXBITS, Esize,Fsize,Asize, NODEBITS );
   }

   c= calloc( N/nnodes*COMPLEXBYTES , 1 );             // allocate amplitudes
   BUFFERSIZE= (1<<18); // number of coefficients in buffer, COMPLEXBYTES bytes each 
   if( BUFFERSIZE> N/nnodes ) BUFFERSIZE=N/nnodes;
   buffer= calloc( BUFFERSIZE*COMPLEXBYTES, 1);    // for communication

   qcprogram(); // this is where you type the quantum program

   clock_gettime(CLOCK_REALTIME,&tim1);
   time= 1.0*(tim1.tv_sec-tim0.tv_sec)+1.e-9*(tim1.tv_nsec-tim0.tv_nsec);
   if(inode==0) printf("\ntime(s):%g\n",time);
   MPI_Finalize();
   exit(0);
}
///////////////////////////////////////////////////////////////////////////////////
