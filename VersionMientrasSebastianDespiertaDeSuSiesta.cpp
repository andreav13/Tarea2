#include <iostream>
#include <fstream> 
#include <cmath> 
#include "Vector.h"

using namespace std;

const double G=1;
const int N=2;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
class Cuerpo;
class Colisionador;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R;
  
public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};

  friend class Colisionador;
};


//Clase Colisionador
class Colisionador{
private:

public:
  void CalculeTodasLasFuerzas(Cuerpo* Planeta);
  void CalculeLaFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
};
  
//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0; R=R0;
}


void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}


void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}


void Cuerpo::Mueva_r(double dt, double Constante){
  r+=V*(Constante*dt);
}


void Cuerpo::Mueva_V(double dt, double Constante){
  V+=F*(Constante*dt)/m;
}

//Funciones de la clase Colisionador
void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Planeta){
  int i,j;
  for(i=0;i<N;i++){
    Planeta[i].BorreFuerza();
  }
  //Calcular todas las fuerzas entre parejas de planetas
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      CalculeLaFuerzaEntre(Planeta[i], Planeta[j]);
    }
  }
}


void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D F1, dr = Planeta2.r-Planeta1.r;
  double aux = G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  F1 = aux*dr;
  Planeta1.AgregueFuerza(F1); Planeta2.AgregueFuerza(F1*(-1));
}



int main(void){
  
  int i;
  double t, dt=50;
  Cuerpo Planeta[N];
  Colisionador Newton;

  double m0=1047, m1=1, r=1000, R0=60, R1=10;
  double M=m0+m1;
  double x0=-m1*r/M, x1=x0+r;

  double omega, Vy0, Vy1, T, tmax;
  omega=sqrt(G*M*pow(r,-3)); Vy0=omega*x0; Vy1=omega*x1; T=2*M_PI/omega; tmax=20*T;

    
  //            (x0, y0, z0, Vx0, Vy0, Vz0, m0, R0)
  Planeta[0].Inicie(x0, 0, 0, 0, Vy0, 0, m0, R0);
  Planeta[1].Inicie(x1, 0, 0, 0, Vy1, 0, m1, R1);

  for (t=0;t<tmax;t+=dt){
    
    cout<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<" "<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;
  
    //Muevase con Omelyan FR.
    for(i=0;i<N;i++){
      Planeta[i].Mueva_r(dt, Zeta);
    }
    Newton.CalculeTodasLasFuerzas(Planeta);
    
    for(i=0;i<N;i++){
      Planeta[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      Planeta[i].Mueva_r(dt, Xi);
    }
   
    Newton.CalculeTodasLasFuerzas(Planeta);
    
    for(i=0;i<N;i++){
      Planeta[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      Planeta[i].Mueva_r(dt, 1-2*(Xi+Zeta)); 
    }
    Newton.CalculeTodasLasFuerzas(Planeta);
    
    for(i=0;i<N;i++){
      Planeta[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      Planeta[i].Mueva_r(dt, Xi);
    }
    
    Newton.CalculeTodasLasFuerzas(Planeta);
    
    for(i=0;i<N;i++){
      Planeta[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      Planeta[i].Mueva_r(dt, Zeta);
    }

    
  }
  
  return 0;

}
