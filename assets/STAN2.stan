// The functions are based on the VineCopula package found on https://github.com/tnagler/VineCopula

functions{
  
  real quantile_t(real u, real nu) {
    real y;
    real sign;
    
    if(u < 0.5){
      y=2 * u;
      sign=-1;
    }else{
      y=2 * (1 - u);
      sign=1;
    }
    
    if(u==0.5){
      return 0;
    }else{
      return sign * sqrt(nu / inv_inc_beta(nu / 2, 0.5, y) - nu);
    }
  }
  
  real dbb1(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t16; real t17; real t38; real t39; real t4; real t5; real t6; real t7; real t9; real t10; real t12; real t13; real t20; real t24; real t25; real t27; real t29; real t32; real t33; real t34; real t36; real t43; real t59;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=pow(u,(-th));
    t2=t1-1.0;
    t3=pow(t2,de);
    t16=1./u;
    t17=1./t2;
    t38=t1*t16;
    t39=t38*t17;
    t4=pow(v,(-th));
    t5=t4-1.0;
    t6=pow(t5,de);
    t7=t3+t6;
    t9=pow(t7,(1./de));
    t10=1.0+t9;
    t12=pow(t10,(-1./th));
    t13=t12*t9;
    t20=1./t10;
    t24=t9*t9;
    t25=t12*t24;
    t27=1./v;
    t29=1./t5;
    t32=t7*t7;
    t33=1./t32;
    t34=t10*t10;
    t36=t33/t34;
    t43=t4*th;
    t59=t43*t27*t29;
    
    out=t25*t6*t27*t4*t29*t36*t3*t39-t13*t6*t43*t27*t29*t33*t3*t38*t17*t20+
    t13*t3*t38*t17*t33*t20*t6*de*t59+t25*t3*t39*t36*t6*t59;
    
    return out;
  }
  
  
  real dbb6(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t4; real t5; real t12; real t16; real t32;
    real t38; real t39; real t40; real t47; real t50; real t61; real t90;
    real t6; real t7; real t8; real t9; real t10; real t11; real t13; real t14;
    real t35; real t36; real t37; real t42; real t48; real t53; real t56;
    real t57; real t59; real t78; real t80; real t87; real t93;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=1.0-u;
    t2=pow(t1,th);
    t3=1.0-t2;
    t4=log(t3);
    t5=pow(-t4,de);
    t12=1/de;
    t16=1/th;
    t32=de-1.0;
    t38=2.0*de;
    t39=-1.0+t38;
    t40=pow(-t4,t39);
    t47=3.0*de-1.0;
    t50=pow(-t4,t32);
    t61=pow(-t4,t47);
    t90=pow(-t4,t38);
    t6=1.0-v;
    t7=pow(t6,th);
    t8=1.0-t7;
    t9=log(t8);
    t10=pow(-t9,de);
    t11=t5+t10;
    t13=pow(t11,t12);
    t14=exp(-t13);
    t35=pow(t11,-2.0*t32*t12);
    t36=t35*th;
    t37=exp(t13);
    t42=pow(-t9,t39);
    t48=pow(-t9,t47);
    t53=t13*de;
    t56=pow(-t9,t32);
    t57=t37*t50*t56;
    t59=t13*th;
    t78=t37-1.0;
    t80=pow(t78*t14,t16);
    t87=t78*t78;
    t93=pow(-t9,t38);
    
    out=(2.0*t36*t37*t40*t42+t36*t37*t48*t50+t53*th*t57-t59*t57+
    t36*t37*t61*t56-2.0*t35*t40*t42-t35*t61*t56-t53*th*t50*t56+t59*t50*t56-
    t35*t48*t50) *t80*t7*t2/t3/t8/t87/(t90+2.0*t5*t10+t93)/t1/t6;
    
    return out;
  }
  
  
  real dbb7(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t4; real t5; real t6; real t7; real t8; real t9; real t11; real t12; real t14; real t15; real t16; real t18; real t20; real t23; real t24; real t25; real t27; real t30; real t31; real t32; real t35; real t37; real t42; real t54;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=1.0-u;
    t2=pow(t1,th);
    t3=1.0-t2;
    t4=pow(t3,-de);
    t5=1.0-v;
    t6=pow(t5,th);
    t7=1.0-t6;
    t8=pow(t7,-de);
    t9=t4+t8-1.0;
    t11=pow(t9,-1.0/de);
    t12=1.0-t11;
    t14=pow(t12,1.0/th);
    t15=t11*t11;
    t16=t14*t15;
    t18=1./t5;
    t20=1./t7;
    t23=t9*t9;
    t24=1./t23;
    t25=t12*t12;
    t27=t24/t25;
    t30=t2/t1;
    t31=1./t3;
    t32=t30*t31;
    t35=t14*t11;
    t37=t6*th;
    t42=1./t12;
    t54=t37*t18*t20;
    
    out=-t16*t8*t6*t18*t20*t27*t4*t32 + t35*t8*t37*t18*t20*t24*t4*t30*t31*t42+
    t35*t4*t30*t31*t24*t42*t8*de*t54+t16*t4*t32*t27*t8*t54;
    
    return out;
  }
  
  
  real dbb8(real u, real v, vector param)
  {
    real th; real de;
    real t2; real t3; real t12; real t16; real t6; real t7; real t10; real t11; real t33; real t38; real t39; real t49; real t59; real t69; real t25; real t26; real t29; real t44; real t45; real t50; real t54; real t62; real t67;
    real out;
    
    th=param[1];
    de=param[2];
    
    t2=1.0-de*u;
    t3=pow(t2,th);
    t10=1.0-de;
    t11=pow(t10,th);
    t12=1.0-t11;
    t16=1/th;
    t33=th*t3;
    t38=2.0*th;
    t39=pow(t10,t38);
    t49=pow(t2,t38);
    t59=pow(t10,3.0*th);
    t69=t12*t12;
    t6=1.0-de*v;
    t7=pow(t6,th);
    t25=t3*t7;
    t26=t11-t7-t3+t25;
    t29=pow(-t26/t12,t16);
    t44=pow(t6,t38);
    t45=t3*t44;
    t50=t49*t7;
    t54=t49*t44;
    t62=-2.0*t25*t11+t25-t33*t7+3.0*t33*t7*t11-3.0*t33*t7*t39+t25*t39+
    2.0* t45*t11-t45*t39+2.0*t50*t11-t50*t39-2.0*t54*t11+t54*t39+t54-
    t50-t45+t33*t7*t59;
    t67=t26*t26;
    out=-de*t29*t62/t6/t2/t67/t69;
    
    return out;
  }
  
  real ta2(real t, real par, real par2, real par3)	//for PDF
  {
    real out;
    real t1; real t2;
    
    t1=pow(par3*t,par);
    t2=pow(par2*(1.0-t),par);
    out=t1+t2;
    return out;
  }
  
  real Tawn2(real t, real par, real par2, real par3)		//for PDF
  {
    real out;
    real t1; real t2; real t3; real t4;
    t1=(1.0-par2)*(1.0-t);
    t2=(1.0-par3)*t;
    t3=ta2(t, par, par2, par3);
    t4=pow(t3,1.0/(par));
    out=t1+t2+t4;
    return out;
  }
  
  real TawnC(real u, real v, real par, real par2, real par3)	// CDF for PDF
  {
    real out;
    real w;
    real A;
    w=log(v)/log(u*v);
    A=Tawn2(w, par, par2, par3);
    out=pow(u*v,A);
    return out;
  }
  
  real d1ta(real t, real par, real par2, real par3)	//for PDF
  {
    real out;
    real t1; real t2;
    t1=par3 * pow((par3*t),par-1.0);
    t2=par2 * pow(par2*(1.0-t),par-1.0);
    out=par*(t1-t2);
    return out;
  }
  
  real d2ta(real t, real par, real par2, real par3)	//for PDF
  {
    real out;
    real t1; real t2;
    t1=pow(par3,2) * pow(par3*t,par-2.0);
    t2=pow(par2,2) * pow(par2*(1.0-t),par-2.0);
    out=par*(par-1) * (t1 + t2);
    return out;
  }
  
  real d1Tawn(real t, real par, real par2, real par3)	//for PDF
  {
    real out;
    real t2; real t1;
    t1=ta2(t, par, par2, par3);
    t2=d1ta(t, par, par2, par3);
    out=par2-par3+1.0/(par) * pow(t1,(1.0/(par)-1.0)) * t2;
    return out;
  }
  
  real d2Tawn(real t, real par, real par2, real par3)	//for PDF
  {
    real out;
    real t2; real t1; real t3;
    
    t1=ta2(t, par, par2, par3);
    t2=d1ta(t, par, par2, par3);
    t3=d2ta(t, par, par2, par3);
    out=1.0/(par) * ( (1.0/(par)-1.0) * pow(t1,(1.0/(par)-2)) * pow(t2,2) + pow(t1,(1.0/(par)-1)) * t3);
    return out;
  }
  
  real dA_du(real u, real v, real par, real par2, real par3)
  {
    real out;
    real dA; real dw; real w;
    w=log(v) / log(u*v);
    dw=-log(v) / (u*pow(log(u*v),2.0));
    dA=d1Tawn(w, par, par2, par3);
    out=dA*dw;
    return out;
  }
  
  real dA_dv(real u, real v, real par, real par2, real par3)
  {
    real out;
    real dA; real dw; real w;
    w=log(v)/log(u*v);
    dw=1.0 / (v*log(u*v)) - log(v) / (v*pow(log(u*v),2));
    dA=d1Tawn(w, par, par2, par3);
    out=dA*dw;
    return out;
  }
  
  real dC_du(real u, real v, real par, real par2, real par3)
  {
    real out;
    real w; real A; real C; real dA;
    w=log(v)/log(u*v);
    A=Tawn2(w, par, par2, par3);
    C=TawnC(u, v, par, par2, par3);
    dA=dA_du(u, v, par, par2, par3);
    out=C * (1.0 / u * A + log(u*v)*dA);
    return out;
  }
  
  real dA_dudv(real u, real v, real par, real par2, real par3)
  {
    real out;
    real dA; real dw_dv; real dw_du; real w; real d2w_dudv; real d2A;
    w=log(v)/log(u*v);
    dw_du=-log(v)/(u*pow(log(u*v),2));
    dw_dv=1.0 / (v*log(u*v)) - log(v) / (v*pow(log(u*v),2));
    d2w_dudv=2*log(v) / (v*u*pow(log(u*v),3)) - 1.0 / (v*u*pow(log(u*v),2));
    dA=d1Tawn(w, par, par2, par3);
    d2A=d2Tawn(w, par, par2, par3);
    out=d2A*dw_dv*dw_du + dA*d2w_dudv;
    return out;
  }
  
  real TawnPDF(real u, real v, real par, real par2, real par3)
  {
    real out;
    real w; real A; real dC; real t3; real t4; real t1; real C; real t5; real t2;
    w=log(v)/log(u*v);
    A=Tawn2(w, par, par2, par3);
    dC=dC_du(u, v, par, par2, par3);
    t3=dA_du(u, v, par, par2, par3);
    t4=dA_dv(u, v, par, par2, par3);
    t1=dC * (1.0/v * A + log(u*v) * t4);
    C=TawnC(u, v, par, par2, par3);
    t5=dA_dudv(u, v, par, par2, par3);
    t2=C * (1.0/u*t4 + 1.0/v*t3 + log(u*v)*t5);
    out=t1+t2;
    return out;
  }
  
  //////////////////////////////////////////////////////////////
  // Function to compute log-likelihood for bivariate copula
  // Input:
  // family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
  // u         first variable of data set
  // v         second variable of data set
  // theta     dependency parameter
  // nu        degrees-of-freedom for students copula
  // loglik    output
  //////////////////////////////////////////////////////////////
  real LL(int family, real u_old, real v_old, real theta, real nu)
  {
    real UMAX=1-1e-12;
    real UMIN=1e-12;
    real XEPS=1e-4;
    real XINFMAX=1.79769e+308;
    real DBLMIN=2.22507e-308;
    vector[2] dat; real rho; real ll=0.0; real t1=0.0; real t2=0.0; real f; vector[2] temp;
    int k=1;
    
    real u=u_old;
    real v=v_old;
    
    if(u<UMIN){
      u=UMIN;
    }else if(u>UMAX){
      u=UMAX;
    }
    if(v<UMIN){
      v=UMIN;
    }else if(v>UMAX){
      v=UMAX;
    }
    
    //Compute log-likelihood:
    if(family==0){
      //independent
      ll=0;
    }else if(family==1){
      //Gaussian
      rho=theta;
      dat[1]=u; dat[2]=v;
      t1=inv_Phi(dat[1]); t2=inv_Phi(dat[2]);
      f=1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==2){
      //Student
      rho=theta;
      dat[1]=u; dat[2]=v;
      t1=quantile_t(dat[1],nu); t2=quantile_t(dat[2],nu);
      f=1/(2*pi()*sqrt(1.0-pow(rho,2.0))*exp(student_t_lpdf(t1|nu,0,1))*exp(student_t_lpdf(t2|nu,0,1)))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==3){
      //Clayton
      if(theta==0){
        ll=0;
      }else if(theta < 1e-10){
        ll=0;
      }else{
        dat[1]=u; dat[2]=v;
        f=log1p(theta)-(1.0+theta)*log(dat[1]*dat[2])-(2.0+1.0/(theta))*log(pow(dat[1],-theta)+pow(dat[2],-theta)-1.0);
        if(f>XINFMAX){
          ll=log(XINFMAX);
        }else if(f<log(DBLMIN)){
          ll=log(DBLMIN);
        }else{
          ll=f;
        }
      }
    }else if(family==4){
      //Gumbel
      dat[1]=u; dat[2]=v;
      t1=pow(-log(dat[1]),theta)+pow(-log(dat[2]),theta);
      f=-pow(t1,1.0/(theta))+(2.0/(theta)-2.0)*log(t1)+(theta-1.0)*log(log(dat[1])*log(dat[2]))-log(dat[1]*dat[2])+log1p((theta-1.0)*pow(t1,-1.0/(theta)));
      if(f>XINFMAX){
        ll=log(XINFMAX);
      }else if(f<log(DBLMIN)){
        ll=log(DBLMIN);
      }else{
        ll=f;
      }
    }else if(family==5){
      // Frank
      if(abs(theta) < 1e-10) {
        ll=0;
      }else {
        dat[1]=u; dat[2]=v;
        f=(theta*(exp(theta)-1.0)*exp(theta*dat[2]+theta*dat[1]+theta))/pow(exp(theta*dat[2]+theta*dat[1])-exp(theta*dat[2]+theta)-exp(theta*dat[1]+theta)+exp(theta),2.0);
        if(log(f)>XINFMAX){
          ll=log(XINFMAX);
        }else if(f < DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(f);
        }
      }
    }else if(family==6){
      //Joe
      f=pow(pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta),1/(theta)-2)*pow(1-u,theta-1)*pow(1-v,theta-1)*(theta-1+pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta));
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==7){
      //BB1
      if(theta==0){
        dat[1]=u; dat[2]=v;
        t1=pow(-log(dat[1]),nu)+pow(-log(dat[2]),nu);
        f=-pow(t1,1/(nu))+(2/(nu)-2)*log(t1)+(nu-1)*log(log(dat[1])*log(dat[2]))-log(dat[1]*dat[2])+log(1+(nu-1)*pow(t1,-1.0/(nu)));
        if(f>XINFMAX){
          ll=log(XINFMAX);
        }else if(f<log(DBLMIN)){
          ll=log(DBLMIN);
        }else{
          ll=f;
        }
      }else{
        vector[2] param;
        real fuc;
        param[1]=theta;
        param[2]=nu;
        fuc=dbb1(u, v, param);
        if(fuc==positive_infinity() || is_nan(fuc)){
          fuc=1;
        }
        if(log(fuc)>XINFMAX){
          ll=log(XINFMAX);
        }else if(fuc<DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(fuc);
        }
      }
    }else if(family==8){
      //BB6
      vector[2] param;
      real fuc;
      param[1]=theta;
      param[2]=nu;
      fuc=dbb6(u, v, param);
      if(fuc==positive_infinity() || is_nan(fuc)){
        fuc=1;
      }
      
      if(log(fuc)>XINFMAX){
        ll=log(XINFMAX);
      }else if(fuc<DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(fuc);
      }
    }else if(family==9){
      //BB7
      if(nu==0){
        f=pow(pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta),1/(theta)-2)*pow(1-u,theta-1)*pow(1-v,theta-1)*(theta-1+pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta));
        if(log(f)>XINFMAX){
          ll=log(XINFMAX);
        }else if(f < DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(f);
        }
      }else{
        vector[2] param;
        real fuc;
        param[1]=theta;
        param[2]=nu;
        fuc=dbb7(u, v, param);
        if(fuc==positive_infinity() || is_nan(fuc)){
          fuc=1;
        }
        
        if(log(fuc)>XINFMAX){
          ll=log(XINFMAX);
        }else if(fuc<DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(fuc);
        }
      }
    }else if(family==10){
      //BB8
      vector[2] param;
      real fuc;
      param[1]=theta;
      param[2]=nu;
      fuc=dbb8(u, v, param);
      if(fuc==positive_infinity() || is_nan(fuc)){
        fuc=1;
      }
      
      if(log(fuc)>XINFMAX){
        ll=log(XINFMAX);
      }else if(fuc<DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(fuc);
      }
    }else if(family==13){
      //rotated Clayton (180 degree)
      if(theta==0){
        ll=0;
      }else if(theta < XEPS){
        ll=0;
      }else{
        dat[1]=1-u; dat[2]=1-v;
        f=(1.0+theta)*pow(dat[1]*dat[2],-1.0-theta)*pow(pow(dat[1],-theta)+pow(dat[2],-theta)-1.0,-2.0-1.0/(theta));
        temp[1]=f;
        temp[2]=0;
        f=max(temp);
        if(log(f)>XINFMAX){
          ll=log(XINFMAX);
        }else if(f < DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(f);
        }
      }
    }else if(family==14){
      //rotated Gumbel (180 degree)
      dat[1]=1-u; dat[2]=1-v;
      t1=pow(-log(dat[1]),theta)+pow(-log(dat[2]),theta);
      t2=exp(-pow(t1,1.0/(theta)));
      f=t2/(dat[1]*dat[2])*pow(t1,-2.0+2.0/(theta))*pow(log(dat[1])*log(dat[2]),theta-1.0)*(1.0+(theta-1.0)*pow(t1,-1.0/(theta)));
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==16){
      //rotated Joe (180 degree)
      u=1-u; v=1-v;
      f=pow(pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta),1/(theta)-2)*pow(1-u,theta-1)*pow(1-v,theta-1)*(theta-1+pow(1-u,theta)+pow(1-v,theta)-pow(1-u,theta)*pow(1-v,theta));
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
      u=1-u; v=1-v;
    }else if(family==17){
      //rotated BB1
      if(theta==0){
        dat[1]=1-u; dat[2]=1-v;
        t1=pow(-log(dat[1]),nu)+pow(-log(dat[2]),nu);
        f=-pow(t1,1/(nu))+(2/(nu)-2)*log(t1)+(nu-1)*log(log(dat[1])*log(dat[2]))-log(dat[1]*dat[2])+log(1+(nu-1)*pow(t1,-1.0/(nu)));
        if(f>XINFMAX){
          ll=log(XINFMAX);
        }else if(f<log(DBLMIN)){
          ll=log(DBLMIN);
        }else{
          ll=f;
        }
      }else{
        vector[2] param;
        real fuc;
        param[1]=theta;
        param[2]=nu;
        
        fuc=dbb1(1-u, 1-v, param);
        
        if(fuc==positive_infinity() || is_nan(fuc)){
          fuc=1;
        }
        
        if(log(fuc)>XINFMAX){
          ll=log(XINFMAX);
        }else if(fuc<DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(fuc);
        }
      }
    }else if(family==18){
      //rotated BB6
      vector[2] param;
      real fuc;
      param[1]=theta;
      param[2]=nu;
      
      fuc=dbb6(1-u,1-v,param);
      
      if(fuc==positive_infinity() || is_nan(fuc)){
        fuc=1;
      }
      
      if(log(fuc)>XINFMAX){
        ll=log(XINFMAX);
      }else if(fuc<DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(fuc);
      }
    }else if(family==19){
      //rotated BB7
      if(nu==0){
        f=pow(pow(u,theta)+pow(v,theta)-pow(u,theta)*pow(v,theta),1/(theta)-2)*pow(u,theta-1)*pow(v,theta-1)*(theta-1+pow(u,theta)+pow(v,theta)-pow(u,theta)*pow(v,theta));
        if(log(f)>XINFMAX){
          ll=log(XINFMAX);
        }else if(f < DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(f);
        }
      }else{
        vector[2] param;
        real fuc;
        param[1]=theta;
        param[2]=nu;
        
        fuc=dbb7(1-u, 1-v, param);
        
        if(fuc==positive_infinity() || is_nan(fuc)){
          fuc=1;
        }
        
        if(log(fuc)>XINFMAX){
          ll=log(XINFMAX);
        }else if(fuc<DBLMIN){
          ll=log(DBLMIN);
        }else{
          ll=log(fuc);
        }
      }
    }else if(family==20){
      //rotated BB8
      vector[2] param;
      real fuc;
      param[1]=theta;
      param[2]=nu;
      
      fuc=dbb8(1-u, 1-v, param);
      
      if(fuc==positive_infinity() || is_nan(fuc)){
        fuc=1;
      }
      
      if(log(fuc)>XINFMAX){
        ll=log(XINFMAX);
      }else if(fuc<DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(fuc);
      }
    }else if(family==104){
      //New: Tawn
      real par3=1.0;
      f=TawnPDF(u, v, theta, nu, par3);
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==114){
      //New: rotated Tawn
      real par3=1.0;
      dat[1]=1-u; dat[2]=1-v;
      f=TawnPDF(dat[1], dat[2], theta, nu, par3);
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==204){
      //New: Tawn2
      real par2=1.0;
      f=TawnPDF(u, v, theta, par2, nu);
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }else if(family==214){
      //New: rotated Tawn2
      real par2=1.0;
      dat[1]=1-u; dat[2]=1-v;
      f=TawnPDF(dat[1], dat[2], theta, par2, nu);
      if(log(f)>XINFMAX){
        ll=log(XINFMAX);
      }else if(f < DBLMIN){
        ll=log(DBLMIN);
      }else{
        ll=log(f);
      }
    }
    
    return ll;
  }
  
  real LL_mod2(int family, real u_old, real v_old, real theta, real nu)
  {
    real UMAX=1-1e-12;
    real UMIN=1e-12;
    real XEPS=1e-4;
    
    real negv;
    real negu;
    real ntheta; real nnu;
    int nfamily;
    
    real u;
    real v;
    
    real loglik;
    
    u=u_old;
    v=v_old;
    
    ntheta=-theta;
    nnu=-nu;
    
    if(u<UMIN){
      u=UMIN;
    }else if(u>UMAX){
      u=UMAX;
    }
    if(v<UMIN){
      v=UMIN;
    }else if(v>UMAX){
      v=UMAX;
    }
    
    if((family==23) || (family==24) || (family==26) || (family==27) || (family==28) || (family==29) || (family==30)){
      // 90 degree rotated copulas
      nfamily=(family)-20;
      negu=1 - u;
      loglik=LL(nfamily, negu,  v, ntheta, nnu);
    }else if((family==33) || (family==34) || (family==36) || (family==37) || (family==38) || (family==39) || (family==40)){
      // 270 degree rotated copulas
      nfamily=(family)-30;
      negv=1 - v;
      loglik=LL(nfamily, u,  negv, ntheta, nnu);
    }else if((family==124) || (family==224)){
      nfamily=(family)-20;
      negu=1 - u;
      loglik=LL(nfamily, v, negu, ntheta, nu);
    }else if((family==134) || (family==234)){
      nfamily=(family)-30;
      negv=1 - v;
      loglik=LL(nfamily, negv, u, ntheta, nu);
    }else{
      loglik=LL(family, u,  v, theta, nu);
    }
    
    return loglik;
  }
  
  // h-func for BB1
  
  real pcondbb1(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t16; real t17; real t4; real t5; real t6; real t7; real t9; real t10; real t12; real t13; real t20;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=pow(u,-th);
    t2=t1-1.;
    t3=pow(t2,de);
    t16=1./u;
    t17=1./t2;
    t4=pow(v,-th);
    t5=t4-1.;
    t6=pow(t5,de);
    t7=t3+t6;
    t9=pow(t7,1/de);
    t10=1.0+t9;
    t12=pow(t10,-1/th);
    t13=t12*t9;
    t20=1./t10;
    out=t13*t3*t1*t16*t17/t7*t20;
    
    return out;
  }
  
  real pcondbb6(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t4; real t5; real t12; real t16; real t6; real t7; real t8; real t9; real t10; real t11; real t13; real t14; real t15; real t17;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=1.0-u;
    t2=pow(t1,th);
    t3=1.0-t2;
    t4=log(t3);
    t5=pow(-t4,de);
    t12=1/de;
    t16=1/th;
    t6=1.0-v;
    t7=pow(t6,th);
    t8=1.0-t7;
    t9=log(t8);
    t10=pow(-t9,de);
    t11=t5+t10;
    t13=pow(t11,t12);
    t14=exp(-t13);
    t15=1.0-t14;
    t17=pow(t15,t16);
    
    out=-t17*t13*t5*t2/t1/t3/t4/t11*t14/t15;
    
    return out;
  }
  
  real pcondbb7(real u, real v, vector param)
  {
    real th; real de;
    real t1; real t2; real t3; real t4; real t6; real t8; real t9; real t11; real t12; real t14;
    real out;
    
    th=param[1];
    de=param[2];
    
    t1=1.0-u;
    t2=pow(t1,1.0*th);
    t3=1.0-t2;
    t4=pow(t3,-1.0*de);
    t6=pow(1.0-v,1.0*th);
    t8=pow(1.0-t6,-1.0*de);
    t9=t4+t8-1.0;
    t11=pow(t9,-1.0/de);
    t12=1.0-t11;
    t14=pow(t12,1.0/th);
    
    out=t14*t11*t4*t2/t1/t3/t9/t12;
    
    return out;
  }
  
  real pcondbb8(real u, real v, vector param)
  {
    real th; real de;
    real t2; real t3; real t12; real t16; real t6; real t7; real t8; real t10; real t11; real t13; real t15; real t17;
    real out;
    
    th=param[1];
    de=param[2];
    
    t2=1.0-de*u;
    t3=pow(t2,th);
    t10=1.0-de;
    t11=pow(t10,th);
    t12=1.0-t11;
    t13=1/t12;
    t16=1/th;
    t6=1.0-de*v;
    t7=pow(t6,th);
    t8=1.0-t7;
    t15=1.0-(1.0-t3)*t8*t13;
    t17=pow(t15,t16);
    
    out=t17*t3/t2*t8*t13/t15;
    
    return out;
  }
  
  //////////////////////////////////////////////////////////////
  // Function to compute h-function for vine simulation and estimation
  // Input:
  // family   copula family (0=independent,  1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
  // u        variable for which h-function computes conditional distribution function
  // v        variable on which h-function conditions
  // theta    parameter for the copula family
  // nu       degrees-of-freedom for the students copula
  // out      output
  //////////////////////////////////////////////////////////////
  real Hfunc(int family, real u_old, real v_old, real theta, real nu)
  {
    real h;
    vector[2] temp2;
    vector[2] temp1;
    real UMAX=1-1e-12;
    real UMIN=1e-12;
    real XEPS=1e-4;
    real x;
    real out;
    
    real u=u_old;
    real v=v_old;
    
    if((v==0) || (u==0)){
      h=0;
    }else if(v==1){
      h=u;
    }else{
      if(family==0){
        //independent
        h=u;
      }else if(family==1){
        //gaussian
        x=(inv_Phi(u) - theta*inv_Phi(v))/sqrt(1.0-pow(theta,2.0));
        if(x<positive_infinity()){
          h=Phi(x);
        }else if((inv_Phi(u) - theta*inv_Phi(v)) < 0){
          h=0;
        }else{
          h=1;
        }
      }else if(family==2){
        //student
        real t1=0.0; real t2=0.0; real mu; real sigma2;
        t1=quantile_t(u,nu); t2=quantile_t(v,nu);
        mu=theta*t2; sigma2=((nu+t2*t2)*(1.0-theta*(theta)))/(nu+1.0);
        h=student_t_cdf((t1-mu)/sqrt(sigma2)|nu+1.0,0,1);
      }else if(family==3){
        //clayton
        if(theta==0){
          h=u;
        }
        if(theta < XEPS){
          h=u;
        }else{
          x=pow(u,-theta)+pow(v,-theta)-1.0;
          h=pow(v,-theta-1.0)*pow(x,-1.0-1.0/(theta));
          if(theta < 0){
            if(x < 0) h=0;
          }
        }
      }else if(family==4){
        //gumbel
        if(theta==1){
          h=u;
        }else{
          h=-(exp(-pow(pow(-log(v),theta)+pow(-log(u),theta),1.0/(theta)))*pow(pow(-log(v),theta)+pow(-log(u),theta),1.0/(theta)-1.0)*pow(-log(v),theta))/(v*log(v));
        }
      }else if(family==5){
        //frank
        if(theta==0){
          h=u;
        }else{
          h=-(exp(theta)*(exp(theta*u)-1.0))/(exp(theta*v+theta*u)-exp(theta*v+theta)-exp(theta*u+theta)+exp(theta));
        }
      }else if(family==6){
        //joe
        if(theta==1){
          h=u;
        }else{
          h=pow(pow(1.0-u,theta) + pow(1.0-v,theta) - pow(1.0-u,theta)*pow(1.0-v,theta),1.0/(theta)-1) * pow(1.0-v,theta-1.0)*(1-pow(1-u,theta));
        }
      }else if(family==7){
        //BB1
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(nu==1){
          if(theta==0){
            h=u;
          }else{
            h=pow(pow(u,-theta)+pow(v,-theta)-1,-1/(theta)-1)*pow(v,-theta-1);
          }
        }else if(theta==0){
          h=-(exp(-pow(pow(-log(v),nu)+pow(-log(u),nu),1.0/(nu)))*pow(pow(-log(v),nu)+pow(-log(u),nu),1.0/(nu)-1.0)*pow(-log(v),nu))/(v*log(v));
        }else{
          h=pcondbb1(v,u,param);
        }
      }else if(family==8){
        //BB6
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        if(theta==1){
          if(nu==1){h=u;
          }else{
            h=-(exp(-pow(pow(-log(v),nu)+pow(-log(u),nu),1.0/(nu)))*pow(pow(-log(v),nu)+pow(-log(u),nu),1.0/(nu)-1.0)*pow(-log(v),nu))/(v*log(v));
          }
        }else if(nu==1){
          h=pow(pow(1.0-u,theta) + pow(1.0-v,theta) - pow(1.0-u,theta)*pow(1.0-v,theta),1.0/(theta)-1) * pow(1.0-v,theta-1.0)*(1-pow(1-u,theta));
        }else{
          h=pcondbb6(v,u,param);
        }
      }else if(family==9){
        //BB7
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(theta==1){
          if(nu==0){
            h=u;
          }else{
            h=pow(pow(u,-nu)+pow(v,-nu)-1,-1/(nu)-1)*pow(v,-nu-1);
          }
        }else if(nu==0){
          h=pow(pow(1.0-u,theta) + pow(1.0-v,theta) - pow(1.0-u,theta)*pow(1.0-v,theta),1.0/(theta)-1) * pow(1.0-v,theta-1.0)*(1-pow(1-u,theta));
        }else{
          h=pcondbb7(v,u,param);
        }
      }else if(family==10){
        //BB8
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(nu==0){
          h=u;
        }else if(nu==1){
          if(theta==1){
            h=u;
          }else{
            h=pow(pow(1.0-u,theta) + pow(1.0-v,theta) - pow(1.0-u,theta)*pow(1.0-v,theta),1.0/(theta)-1) * pow(1.0-v,theta-1.0)*(1-pow(1-u,theta));
          }
        }else{
          h=pcondbb8(v,u,param);
        }
      }else if(family==13){
        //rotated clayton (180 degree)
        if(theta==0){
          h=u;
        }if(theta < XEPS){
          h=u;
        }else{
          u=1-u;
          v=1-v;
          x=pow(u,-theta)+pow(v,-theta)-1.0;
          h=pow(v,-theta-1.0)*pow(x,-1.0-1.0/(theta));
          h=1-h;
          u=1-u;
          v=1-v;
        }
      }else if(family==14){
        //rotated gumbel (180 degree)
        v=1-v;
        u=1-u;
        h=-(exp(-pow(pow(-log(v),theta)+pow(-log(u),theta),1.0/(theta)))*pow(pow(-log(v),theta)+pow(-log(u),theta),1.0/(theta)-1.0)*pow(-log(v),theta))/(v*log(v));
        h=1-h;
        u=1-u;
        v=1-v;
      }else if(family==16){
        v=1-v;
        u=1-u;
        h=pow(pow(1.0-u,theta) + pow(1.0-v,theta) - pow(1.0-u,theta)*pow(1.0-v,theta),1.0/(theta)-1) * pow(1.0-v,theta-1.0)*(1-pow(1-u,theta));
        h=1-h;
        u=1-u;
        v=1-v;
      }else if(family==17){
        //rotated BB1
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(nu==1){
          if(theta==0){
            h=u;
          }else{
            h=pow(pow(1-u,-theta)+pow(1-v,-theta)-1,-1/(theta)-1)*pow(1-v,-theta-1);
            h=1-h;
          }
        }else if(theta==0){
          h=-(exp(-pow(pow(-log(1-v),nu)+pow(-log(1-u),nu),1.0/(nu)))*pow(pow(-log(1-v),nu)+pow(-log(1-u),nu),1.0/(nu)-1.0)*pow(-log(1-v),nu))/((1-v)*log(1-v));
          h=1-h;
        }else{
          v=1-v;
          u=1-u;
          h=pcondbb1(v,u,param);
          u=1-u;
          v=1-v;
          h=1-h;
        }
      }
      else if(family==18){
        //rotated BB6
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(theta==1){
          if(nu==1){
            h=u;
          }else{
            h=-(exp(-pow(pow(-log(1-v),nu)+pow(-log(1-u),nu),1.0/(nu)))*pow(pow(-log(1-v),nu)+pow(-log(1-u),nu),1.0/(nu)-1.0)*pow(-log(1-v),nu))/((1-v)*log(1-v));
            h=1-h;
          }
        }else if(nu==1){
          h=pow(pow(u,theta) + pow(v,theta) - pow(u,theta)*pow(v,theta),1.0/(theta)-1) * pow(v,theta-1.0)*(1-pow(u,theta));
          h=1-h;
        }else{
          v=1-v;
          u=1-u;
          h=pcondbb6(v,u,param);
          u=1-u;
          v=1-v;
          h=1-h;
        }
      }else if(family==19){
        //rotated BB7
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(theta==1){
          if(nu==0){
            h=u;
          }else{
            h=pow(pow(1-u,-nu)+pow(1-v,-nu)-1,-1/(nu)-1)*pow(1-v,-nu-1);
            h=1-h;
          }
        }else if(nu==0){
          h=pow(pow(u,theta) + pow(v,theta) - pow(u,theta)*pow(v,theta),1.0/(theta)-1) * pow(v,theta-1.0)*(1-pow(u,theta));
          h=1-h;
        }else{
          v=1-v;
          u=1-u;
          h=pcondbb7(v,u,param);
          u=1-u;
          v=1-v;
          h=1-h;
        }
      }else if(family==20){
        //rotated BB8
        vector[2] param;
        param[1]=theta;
        param[2]=nu;
        
        if(nu==0){
          h=u;
        }else if(nu==1){
          if(theta==1){
            h=u;
          }else{
            h=pow(pow(u,theta) + pow(v,theta) - pow(u,theta)*pow(v,theta),1.0/(theta)-1) * pow(v,theta-1.0)*(1-pow(u,theta));
            h=1-h;
          }
        }else{
          v=1-v;
          u=1-u;
          h=pcondbb8(v,u,param);
          u=1-u;
          v=1-v;
          h=1-h;
        }
      }
    }
    temp1[1]=h;
    temp1[2]=UMAX;
    temp2[1]=min(temp1);
    temp2[2]=UMIN;
    out=max(temp2);
    
    return out;
  }
  
  // Since the h function is not symmetric in case of real Gumbel and real Clayton we have two implement both separately,
  // i.e. Hfunc1 and Hfunc2
  real  Hfunc1(int family, real u_old, real v_old,real theta,real nu)
  {
    real UMAX=1-1e-12;
    real UMIN=1e-12;
    real XEPS=1e-4;
    vector[2] temp1;
    vector[2] temp2;
    real u=u_old;
    real v=v_old;
    real out;
    real negv; real negu;
    real ntheta; real nnu;
    int nfamily;
    
    ntheta=-theta;
    nnu=-nu;
    
    if(u<UMIN){
      u=UMIN;
    }else if(u>UMAX){
      u=UMAX;
    }
    if(v<UMIN){
      v=UMIN;
    }else if(v>UMAX){
      v=UMAX;
    }
    
    if(((family==23) || (family==24) || (family==26) || (family==27) || (family==28) || (family==29) || (family==30))){
      nfamily=(family)-20;
      negv=1 - v;
      out=Hfunc(nfamily, u, negv, ntheta, nnu);
    }else if(((family==33) || (family==34) || (family==36) || (family==37) || (family==38) || (family==39) || (family==40))){
      nfamily=(family)-30;
      negu=1 - u;
      out=Hfunc(nfamily, negu, v, ntheta, nnu);
      out=1-out;
    }else if(family==104){
      real par3=1;
      out=dC_du(v,u,theta,nu,par3);
    }else if(family==114){
      real par3=1;
      negv=1-v;
      negu=1-u;
      out=dC_du(negv,negu,theta,nu,par3);
      out=1-out;
    }else if(family==124){
      real par3=nu;
      real par2=1;
      negv=1-v;
      out=dC_du(negv,u,ntheta,par2,par3);
    }else if(family==134){
      real par3=nu;
      real par2=1;
      negu=1-u;
      out=dC_du(v,negu,ntheta,par2,par3);
      out=1-out;
    }else if(family==204){
      real par3=nu;
      real par2=1;
      out=dC_du(v,u,theta,par2,par3);
    }else if(family==214){
      real par3=nu;
      real par2=1;
      negv=1-v;
      negu=1-u;
      out=dC_du(negv,negu,theta,par2,par3);
      out=1-out;
    }else if(family==224){
      real par3=1;
      negv=1-v;
      out=dC_du(negv,u,ntheta,nu,par3);
    }else if(family==234){
      real par3=1;
      negu=1-u;
      out=dC_du(v,negu,ntheta,nu,par3);
      out=1-out;
    }else{
      out=Hfunc(family, u, v, theta, nu);
    }
    // ensure that results are in [0,1]
    temp1[1]=UMIN;
    temp1[2]=out;
    temp2[1]=UMAX;
    temp2[2]=max(temp1);
    out=min(temp2);
    
    return out;
  }
  
  real Hfunc2(int family, real v_old, real u_old, real theta, real nu)
  {
    real negv; real negu;
    
    real UMAX=1-1e-12;
    real UMIN=1e-12;
    real XEPS=1e-4;
    
    vector[2] temp1;
    vector[2] temp2;
    real u=u_old;
    real v=v_old;
    real out;
    real ntheta; real nnu;
    int nfamily;
    
    ntheta=-theta;
    nnu=-nu;
    
    if(u<UMIN){
      u=UMIN;
    }else if(u>UMAX){
      u=UMAX;
    }
    if(v<UMIN){
      v=UMIN;
    }else if(v>UMAX){
      v=UMAX;
    }
    
    if((family==23) || (family==24) || (family==26) || (family==27) || (family==28) || (family==29) || (family==30)){
      nfamily=family-20;
      negv=1 - v;
      out=Hfunc(nfamily, negv, u, ntheta, nnu);
      out=1-out;
    }else if((family==33) || (family==34) || (family==36) || (family==37) || (family==38) || (family==39) || (family==40)){
      nfamily=(family)-30;
      negu=1 - u;
      out=Hfunc(nfamily, v, negu, ntheta, nnu);
    }else if((family==104) || (family==204) || (family==114) || (family==214)){
      // switch u and v and change type
      if((family)/100==1){
        nfamily=(family) + 100;
      }
      if((family)/100==2){
        nfamily=(family) - 100;
      }
      negu=1 - u;
      out=Hfunc1(nfamily, v, u, theta, nu);
    }else if((family==124) || (family==224) || (family==134) || (family==234)){
      // switch u and v and change type
      if((family)/100==1){
        nfamily=(family) + 100;
      }
      if((family)/100==2){
        nfamily=(family) - 100;
      }
      negv=1 - v;
      negu=1 - u;
      out=Hfunc1(nfamily, negv, negu, theta, nu);
      out=1 - out;
    }else{
      // switch u and v
      out=Hfunc(family, v, u, theta, nu);
    }
    
    // ensure that results are in [0,1]
    temp1[1]=UMIN;
    temp1[2]=out;
    temp2[1]=UMAX;
    temp2[2]=max(temp1);
    out=min(temp2);
    
    return out;
  }
  
  real VineLogLikRvine2(int T, int d, array[] int family, array[] int maxmat, array[] int matri, array[] int condirect,
  array[] int conindirect, vector par, vector par2, vector daten)
  {
    int l; int m;
    
    vector[T] out;
    
    array[d,d] real value2;
    
    array[d,T] real x;
    array[d,d] real vdirect;
    array[d,d] real vindirect;
    array[d,d] real theta;
    array[d,d] real nu;
    array[d,d] int fam;
    
    int i;
    int k;
    
    vector[T] sumsplitlog;
    
    for (t in 1:T)
    {
      sumsplitlog[t]=0;
    }
    
    
    //Initialize
    l=1;
    for(ifor in 1:d)
    {
      for (t in 1:T)
      {
        x[ifor][t]=daten[l];
        l=l+1;
      }
    }
    
    for(ifor in 1:d)
    {
      for(j in 1:d)
      {
        theta[ifor][j]=par[ifor+(d)*(j-1)];
        nu[ifor][j]=par2[ifor+(d)*(j-1)]   ;
        fam[ifor][j]=family[ifor+(d)*(j-1)];
      }
    }
    
    for(t in 1:T)
    {
      for(ifor in 1:d)
      {
        vdirect[d][ifor]=x[d-ifor+1][t];
      }
      
      for(ifor in 1:d-1)
      {
        i=d-ifor;
        for(kfor in (i+1):d)
        {
          k=d+i+1-kfor;
          m=maxmat[k+d*(i-1)];
          if(m==matri[k+d*(i-1)]){
            value2[k][i]=LL_mod2(fam[k][i],vdirect[k][d-m+1],vdirect[k][i],theta[k][i],nu[k][i]);
            
            if(condirect[k+d*(i-1)-1]==1){
              vdirect[k-1][i]=Hfunc1(fam[k][i],vdirect[k][i],vdirect[k][d-m+1],theta[k][i],nu[k][i]);
            }
            if(conindirect[k+d*(i-1)-1]==1){
              vindirect[k-1][i]=Hfunc2(fam[k][i],vdirect[k][(d-m+1)],vdirect[k][i],theta[k][i],nu[k][i]);
            }
          }else{
            value2[k][i]=LL_mod2(fam[k][i],vindirect[k][(d-m+1)],vdirect[k][i],theta[k][i],nu[k][i]);
            
            if(condirect[k+d*(i-1)-1]==1){
              vdirect[k-1][i]=Hfunc1(fam[k][i],vdirect[k][i],vindirect[k][(d-m+1)],theta[k][i],nu[k][i]);
            }
            if(conindirect[k+d*(i-1)-1]==1){
              vindirect[k-1][i]=Hfunc2(fam[k][i],vindirect[k][(d-m+1)],vdirect[k][i],theta[k][i],nu[k][i]);
            }
          }
          
          sumsplitlog[t] += value2[k][i];
        }
      }
      out[t]=sumsplitlog[t];
    }
    
    return sum(out);
    
  }
  
}

///////////
// The data that is prepared in in the document R and given to STAN
///////////

data{
  int T; //T=1
  int dataflip;
  int d;
  array[d] int o;
  int d1;
  int d2;
  array[d*d] int family;
  array[d*d] int maxmat;
  array[d*d] int matri;
  array[d*d] int condirect;
  array[d*d] int conindirect;
  vector[d*d] par;
  vector[d*d] par2;
  vector[d] indexcon; //the indices
  vector<lower=0, upper=1>[d] ucon;
}

///////////
// What we want to sample from
///////////

parameters{
  vector<lower=0, upper=1>[d1] ucalculate;
}

///////////
// The model
///////////

model{
  vector[d] daten;
  vector[d] datennew;
  vector[d] check;
  
  vector[1] temp1;
  vector[1] temp2;
  vector[1] temp3;
  
  int k;
  
  for(i in 1:d)//search for the index i in indexcon
  {
    for(j in 1:d)
    {
      if(indexcon[j]==i)
      {
        daten[i]=ucon[j];
        check[i]=1;
      }
    }
  }
  
  k=1;
  for (i in 1:d)
  {
    if(check[i]!=1)
    {
      daten[i]=ucalculate[k];
      k=k+1;
    }
    else{
      k=k;
    }
  }
  if(dataflip==1){
    for(i in 1:d){
      datennew[i]=daten[o[d+1-i]];
    }
  }
  else{
    datennew=daten;
  }
  
  target += VineLogLikRvine2(T, d, family, maxmat, matri, condirect, conindirect, par, par2, datennew);
  
}
