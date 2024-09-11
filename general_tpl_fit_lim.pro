;fit a general triple-power-law spectrum with energy E, energy error E_err, flux J and flux error J_err
;
;INPUT:
;  E(eV),J(cm^-2*sr^-1*eV^-1),E_err,J_err
;  p0=[lnA,beta1,beta2,beta3,lnE1,lnE2,alpha_1,alpha_2]
;
;EXAMPLE:
;pro fit_tpl_spectrum
;  E1=5d3
;  E2=30d+3
;  beta1=2d
;  beta2=3d
;  beta3=5d
;  A=exp(10d)
;  alpha_1=5d
;  alpha_2=5d
;
;  E=1d2*1d3^(dindgen(18)/17)
;  J=A*E^(-beta1)*(1+(E/E1)^alpha_1)^((beta1-beta2)/alpha_1)*(1+(E/E2)^alpha_2)^((beta2-beta3)/alpha_2)*(1+0.1*randomn(seed,18))
;  E_err=0.2*E
;  J_err=0.1*J
;
;  p0=[20d,2d,5d,8d,7d,11d,5d,5d]
;  general_tpl_fit_lim,E,J,x_err=E_err,y_err=J_err,p=p0,fita=[1,1,1,1,1,1,1,1],tol=tol,chisq=chisq,reduced_chisq=reduced_chisq,Q=Q,flux_fit=flux_fit,itmax=itmax,fita=fita,covar=covar,$
;  a_fit=a_fit,a_std=a_std,beta1_fit=beta1_fit,beta1_std=beta1_std,beta2_fit=beta2_fit,beta2_std=beta2_std,beta3_fit=beta3_fit,beta3_std=beta3_std,$
;  E1_fit=E1_fit,E1_std=E1_std,E2_fit=E2_fit,E2_std=E2_std,alpha1_fit=alpha1_fit,alpha1_std=alpha1_std,alpha2_fit=alpha2_fit,alpha2_std=alpha2_std,$
;  ;  limit=limit,$ ; flag of 5 limits, 1:tpl limit;  2:double power-law; 3:log parabola;  4:exponetial cut off;
;  Ec_fit=Ec_fit,Ec_std=Ec_std,$ ;parameters for exponential cut off
;  gamma1_fit=gamma1_fit,gamma1_std=gamma1_std,gamma2_fit=gamma2_fit,gamma2_std=gamma2_std,$ ;parameters for log parabola approximation
;  w0_fit=w0_fit,w0_std=w0_std,$
;  numerical=numerical,$ ;numerical derivation for d_lny/d_lnx
;  normalization=normalization,$ ;normalization for reduced_chisq
;  ECVI=ECVI
;end

pro general_tpl_fit_lim,x,y,x_err=x_err,y_err=y_err,p=p,tol=tol,chisq=chisq,reduced_chisq=reduced_chisq,Q=Q,flux_fit=flux_fit,itmax=itmax,fita=fita,covar=covar,$
  a_fit=a_fit,a_std=a_std,beta1_fit=beta1_fit,beta1_std=beta1_std,beta2_fit=beta2_fit,beta2_std=beta2_std,beta3_fit=beta3_fit,beta3_std=beta3_std,$
  E1_fit=E1_fit,E1_std=E1_std,E2_fit=E2_fit,E2_std=E2_std,alpha1_fit=alpha1_fit,alpha1_std=alpha1_std,alpha2_fit=alpha2_fit,alpha2_std=alpha2_std,$
  ;  limit=limit,$ ; flag of 4 limits, 1:tpl 2:dpl 3:LP(Log parabola) 4:ER(Exponential cut-off) 5:kappa
  Ec_fit=Ec_fit,Ec_std=Ec_std,$ ;parameters for exponential cut off
  gamma1_fit=gamma1_fit,gamma1_std=gamma1_std,gamma2_fit=gamma2_fit,gamma2_std=gamma2_std,$ ;parameters for log parabola approximation
  w0_fit=w0_fit,w0_std=w0_std,$
  numerical=numerical,$ ;numerical derivation for d_lny/d_lnx
  normalization=normalization,$ ;normalization for reduced_chisq
  ECVI=ECVI
  ;x: energy (eV)
  ;y: flux (cm^-2*sr^-1*eV^-1)
  ;p: [ln(A),beta1,beta2,beta3,ln(E1),ln(E2),alpha1,alpha2]

  p_max=[5d2,1d4,1d5,1d6,1d2,1d2,5d2,5d2]; bundary of p

  if undefined(x_err) then x_err=0.1*x
  if undefined(y_err) then y_err=0.1*y

  ;transfer into log scale
  lnx=alog(x)
  lny=alog(y)
  lnx_err=x_err/x
  lny_err=y_err/y

  if undefined(p) then begin
    ind_sort=sort(x)
    p1_def=-(lny[ind_sort[n_elements(x)/3]]-lny[ind_sort[0]])/(lnx[ind_sort[n_elements(x)/3]]-lnx[ind_sort[0]])
    p2_def=-(lny[ind_sort[2*n_elements(x)/3]]-lny[ind_sort[n_elements(x)/3]])/(lnx[2*n_elements(x)/3]-lnx[ind_sort[n_elements(x)/3]])
    p3_def=-(lny[ind_sort[-1]]-lny[ind_sort[2*n_elements(x)/3]])/(lnx[-1]-lnx[ind_sort[2*n_elements(x)/3]])
    p4_def=lnx[ind_sort[n_elements(x)/3]]
    p5_def=lnx[ind_sort[2*n_elements(x)/3]]
    p6_def=min([2d1/(max(lnx)-min(lnx)),5d])
    p7_def=min([2d1/(max(lnx)-min(lnx)),5d])
    p0_def=lny[0]-general_tpl(lnx[0],p=[0d,p1_def,p2_def,p3_def,p4_def,p5_def,p6_def,p7_def])
    p=[p0_def,p1_def,p2_def,p3_def,p4_def,p5_def,p6_def,p7_def] ;initial guess of p
  endif
  if undefined(itmax) then itmax=2000
  if undefined(tol) then tol=1d-8
  if undefined(fita) then fita=[1,1,1,1,1,1,1,1]

  ;BFGS method
  BFGS:
  g=gt_dchisq_dp(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)
  k=0
  H=identity(8,/double)*(fita#transpose(fita))
  flag=0
  limit=0

  norm_gfita=findgen(itmax)
  while k lt itmax and norm(g*fita) gt tol do begin
    dp=-H#(g*fita)/norm(H#(g*fita))
    alpha=gt_alpha_seek(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,dp=dp,numerical=numerical)
    p=p+alpha*dp
    delta_p=alpha*dp
    delta_g=gt_dchisq_dp(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)$
      -gt_dchisq_dp(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p-alpha*dp,numerical=numerical)
    p_gHg=total(transpose(delta_g)#H#delta_g)
    omiga=p_gHg/2*(delta_p/total(delta_g*delta_p)-H#delta_g/p_gHg)
    H=H-H#delta_g#transpose(delta_g)#H/p_gHg+delta_p#transpose(delta_p)/total(delta_g*delta_p)+omiga#transpose(omiga)
    H=H*(fita#transpose(fita))
    delta_chisq=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)$
      -general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p-alpha*dp,numerical=numerical)

    ;;;;;limit criteria
    if (-delta_chisq/alpha lt 1d-3*norm(g*fita)) and limit eq 0 then begin
      reduced_chisq_gs=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)/(n_elements(x)-total(fita))
      p_tpl=p
      reduced_chisq_tpl=2*reduced_chisq_gs
      p_dpl=p
      reduced_chisq_dpl=2*reduced_chisq_gs
      p_eco=p
      reduced_chisq_eco=2*reduced_chisq_gs
      p_lp=p
      reduced_chisq_lp=2*reduced_chisq_gs
      p_tlp=p
      reduced_chisq_tlp=2*reduced_chisq_gs
      p_lper=p
      reduced_chisq_lper=2*reduced_chisq_gs
      p_kappa=p
      reduced_chisq_kappa=2*reduced_chisq_gs

      ;; tpl limit
      ;; beta1&beta2, beta2&beta3 can be distinguished
      ;; >=2 energy channels between the two transition[E1_rightband,E2_leftband]
      if p[4] gt min(lnx) and p[4] lt max(lnx) and p[5] gt min(lnx) and $
        p[5] lt max(lnx) and abs(p[2]-p[3]) gt 0.2 and abs(p[1]-p[2]) gt 0.2 $
        and n_elements(where(lnx gt (p[4]+2d/p[6]) and lnx lt (p[5]-2d/p[7]))) ge 2 and $
        total(fita) eq 8 then begin
        delta_lnx_1=abs(lnx-p[4])
        delta_lnx_2=abs(lnx-p[5])
        ind_sort_1=sort(delta_lnx_1)
        ind_sort_2=sort(delta_lnx_2)
        p_lower_1=2/delta_lnx_1[ind_sort_1[1]]
        p_lower_2=2/delta_lnx_2[ind_sort_2[1]]
        p_tpl[6]=min([2*p_lower_1,p_max[6]])
        p_tpl[7]= min([2*p_lower_2,p_max[7]])
        reduced_chisq_tpl=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p_tpl,numerical=numerical)/(n_elements(x)-total(fita*[1,1,1,1,1,1,0,0]))
        ;print,'reduced_chisq_tpl='+strmid(strcompress(string(reduced_chisq_tpl),/rem),0,4)
        ;print,'reduced_chisq_gs='+strmid(strcompress(string(reduced_chisq_gs),/rem),0,4)
        ;print,p_tpl
      endif

;      ;; dpl limit
;      ;; beta2&beta1 can't be distinguished or E1 less than Emin
;      ;; beta3&beta2 can't be distinguished or E2 greater than Emax
;      if (abs(p[1]-p[2]) lt 0.1 or abs(p[2]-p[3]) lt 0.1 or p[4] lt 1.1*min(lnx) or p[5] gt max(lnx)) and total(fita) eq 8 then begin
;        if abs(p[1]-p[2]) lt 0.1 or p[4] lt 1.1*min(lnx) then begin
;          p_dpl[4]=0.1*min(lnx)
;          p_dpl[6]=10d
;        endif
;        if abs(p[2]-p[3]) lt 0.1 or p[5] gt max(lnx) then begin
;          p_dpl[5]=p[4]
;          p_dpl[4]=0.1*min(lnx)
;          p_dpl[1]=p[1]
;          p_dpl[2]=p[1]
;          p_dpl[3]=p[2]
;          p_dpl[7]=p[6]
;          p_dpl[6]=10d
;        endif
;        reduced_chisq_dpl=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p_eco,numerical=numerical)/(n_elements(x)-total(fita*[1,0,1,1,0,1,0,1]))
;        print,'reduced_chisq_dpl='+strmid(strcompress(string(reduced_chisq_dpl),/rem),0,4)
;        print,'reduced_chisq_gs='+strmid(strcompress(string(reduced_chisq_gs),/rem),0,4)
;      endif

      ;; logarithmic parabola limit
      ;; alpha1&alpha2 are small enough(0.1)
      if 4/p[6] gt max(lnx)-min(lnx) and 4/p[7] gt max(lnx)-min(lnx) and total(fita) ge 7 then begin
        lnx0_1=min(lnx)*2d/3d + max(lnx)/3d
        lnx0_2=min(lnx)/3d + max(lnx)*2d/3d
        lny0_1=general_tpl(lnx0_1,p=p)
        lny0_2=general_tpl(lnx0_2,p=p)
        lny1=general_tpl(min(lnx),p=p)
        lny2=general_tpl(max(lnx),p=p)
        beta0_1=-(lny0_2-lny1)/(lnx0_2-min(lnx))
        beta0_2=-(lny2-lny0_1)/(max(lnx)-lnx0_1)
        delta_beta1=8/p_lp[6]*(lny1+lny0_2-2*lny0_1)/(lnx0_2-min(lnx))^2
        delta_beta2=8/p_lp[7]*(lny0_1+lny2-2*lny0_2)/(max(lnx)-lnx0_1)^2
        p_lp[6]=0.5/(lnx0_2-min(lnx))          ;(E0/E)^alpha in range [exp(-0.25),exp(0.25)]
        p_lp[7]=0.5/(max(lnx)-lnx0_1)
        print,'beta0_1-delta_beta1='+strmid(strcompress(string(beta0_1-delta_beta1),/rem),0,4)
        print,'beta0_2+delta_beta2='+strmid(strcompress(string(beta0_2+delta_beta2),/rem),0,4)
        p_lp[1]=beta0_1+delta_beta1
        p_lp[2]=-(lny0_2-lny0_1)/(lnx0_2-lnx0_1)
        p_lp[3]=beta0_2-delta_beta2
        p_lp[4]=lnx0_1
        p_lp[5]=lnx0_2
        p_lp[0]=lny0_1+p_lp[1]*lnx0_1+(p_lp[2]-p_lp[1])/p_lp[6]*alog(2)+(p_lp[3]-p_lp[2])/p_lp[7]*alog(1+exp(p_lp[7]*(lnx0_1-lnx0_2)))
        reduced_chisq_lp=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p_lp,numerical=numerical)/(n_elements(x)-total(fita*[1,1,1,1,0,0,0,0]))
        ;print,'reduced_chisq_lp='+strmid(strcompress(string(reduced_chisq_lp),/rem),0,4)
        ;print,'reduced_chisq_gs='+strmid(strcompress(string(reduced_chisq_gs),/rem),0,4)
      endif

      ;; exponential cut-off limit(ER-like limit)
      ;; beta3 is large enough(beta3>beta2,beta3>beta1  weak condition)
      ;; E2 is large enough(E2>Emax  weak condition)
      if p[3] gt p[2] and p[3] gt p[1] and p[5] gt max(lnx) and total(fita) eq 8 then begin
        eco_limit_1=0
        eco_limit_2=0
        Ec_es_2=p[5]+alog(p[7]/(p[3]-p[2]))/p[7]
        p_eco[3]=1d5

        E2_left=p[5]-2d/p[7]
        if (E2_left-p[4]) lt 0.1*(E2_left-min(x)) then begin
          p_eco[1]=p[1]
          p_eco[2]=p[1]
        endif
        if (p[4]-min(x)) lt 0.1*(E2_left-min(x)) then begin
          p_eco[1]=p[2]
          p_eco[2]=p[2]
        endif else begin
          p_eco[1]=(p[1]+p[2])/2d
          p_eco[2]=(p[1]+p[2])/2d
        endelse

        lnx_mid=(min(lnx)+max(lnx))/2d
        lny_x_mid=general_tpl(lnx_mid,p=p)
        lny_E1=general_tpl(p_eco[4],p=p)
        p_new_2=Ec_es_2+alog((p_eco[3]-p_eco[2])/p_eco[7])/p_eco[7]
        p_eco[5]=p_new_2
        p_eco[0]=lny_E1-(-p_eco[1]*p_eco[4]+(p_eco[1]-p_eco[2])/p_eco[6]*alog(1+exp(p_eco[6]*(p_eco[4]-p_eco[4])))+(p_eco[2]-p_eco[3])/p_eco[7]*alog(1+exp(p_eco[7]*(p_eco[4]-p_eco[5]))))
        reduced_chisq_eco=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p_eco,numerical=numerical)/(n_elements(x)-total(fita*[1,1,0,0,0,1,0,1]))
        ;print,'reduced_chisq_eco='+strmid(strcompress(string(reduced_chisq_eco),/rem),0,4)
        ;print,'reduced_chisq_gs='+strmid(strcompress(string(reduced_chisq_gs),/rem),0,4)
      endif

      ;; !!!!
      ;; kappa distribution(beta1=-1)
      ;; beta1<0  and  Emin<E1<Emax
      if p[1] lt 0 and p[4] gt min(lnx) and p[4] lt max(lnx) and total(fita) eq 8 then begin
        p_kappa[3]=30d
        p_kappa[1]=-1d
        p_kappa[6]=1d
        p_kappa[4]=max(lnx)-2d/p_kappa[6] ;;E1
        lny_E1=general_tpl(p[4],p=p)
        lny_Emax=general_tpl(max(lnx),p=p)
        p_kappa[2]=-(lny_Emax-lny_E1)/(max(lnx)-p[4])

        reduced_chisq_kappa=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p_kappa,numerical=numerical)/(n_elements(x)-total(fita*[1,0,1,0,1,1,0,1]))
        ;print,'reduced_chisq_kappa='+strmid(strcompress(string(reduced_chisq_kappa),/rem),0,4)
        ;print,'reduced_chisq_gs='+strmid(strcompress(string(reduced_chisq_gs),/rem),0,4)
      endif

      ;; output limit
      lim1=reduced_chisq_tpl/reduced_chisq_gs
      lim2=reduced_chisq_dpl/reduced_chisq_gs
      lim3=reduced_chisq_lp/reduced_chisq_gs
      lim4=reduced_chisq_eco/reduced_chisq_gs
      lim5=reduced_chisq_kappa/reduced_chisq_gs
      min_l=min([lim1,lim2,lim3,lim4,lim5])
      print,[lim1,lim2,lim3,lim4,lim5]
      if min_l eq lim1 and min_l lt 1.2 then begin
        p=p_tpl
        message,'alpha1(2) greater than '+strmid(strcompress(string(p[6]/2),/rem),0,4)+'('+$
          strmid(strcompress(string(p[7]/2),/rem),0,4)+')'+', triple power-law approximation!',/continue
        limit=1
        fita[6]=0
        fita[7]=0
        H=H*(fita#transpose(fita))
        k=0
      endif

      if min_l eq lim2 and min_l lt 1.2 then begin
        p=p_dpl
        message,'alpha1(2) greater than '+strmid(strcompress(string(p[6]/2),/rem),0,4)+'('+$
          strmid(strcompress(string(p[7]/2),/rem),0,4)+')'+', double power-law approximation!',/continue
        limit=2
        fita[1]=0
        fita[4]=0
        fita[6]=0
        fita[7]=0
        H=H*(fita#transpose(fita))
        k=0
      endif

      if min_l eq lim3 and min_l lt 1.2 then begin
        p=p_lp
        message,'log parabola!',/continue
        limit=3
        fita[2]=0
        fita[4]=0
        fita[5]=0
        fita[6]=0
        fita[7]=0
        H=H*(fita#transpose(fita))
        k=0
      endif

      if min_l eq lim4 and min_l lt 1.2 then begin
        p=p_eco
        message,'exponential cut off!',/continue
        limit=4
        p[4]=0.1*min(lnx)
        p[6]=10d
        fita[1]=0
        fita[3]=0
        fita[4]=0
        fita[6]=0
        H=H*(fita#transpose(fita))
        k=0
      endif

      if min_l eq lim5 and min_l lt 1.2 then begin
        p=p_kappa
        message,'kappa distribution!',/continue
        limit=6
        p[5]=30d
        p[7]=10d
        fita[3]=0
        fita[5]=0
        fita[7]=0
        H=H*(fita#transpose(fita))
        k=0
      endif
    endif

    overflow=0
    for kk=0,7 do begin
      if p[kk] gt p_max[kk] then begin
        overflow=1
        break
      endif
    endfor
    if overflow then begin
      message,'overflow, stop at residual='+strcompress(string(norm(g*fita)),/rem)+'!',/continue
      break
    endif

    if p[4] gt p[5] then begin ;make sure that E2>E1
      p[2]=p[1]+p[3]-p[2]
      temp=[p[4],p[5],p[6],p[7]]
      p[4]=temp[1]
      p[5]=temp[0]
      p[6]=temp[3]
      p[7]=temp[2]
    endif
    g=gt_dchisq_dp(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)
    k++
    if k lt itmax then norm_gfita[k]=norm(g*fita)
  endwhile

;; run for simulation
;  print,'------norm(g*fita)-------'
;  if k lt itmax then begin
;    print,norm_gfita[0:k]
;  endif else begin
;    print,norm_gfita[0:k-1]
;  endelse

  if k eq itmax and norm(g*fita) gt tol then begin
    message,'Fail to converge, stop at residual='+strcompress(string(norm(g*fita)),/rem)+'!',/continue
  endif

  a_fit=p[0]
  beta1_fit=p[1]
  beta2_fit=p[2]
  beta3_fit=p[3]
  E1_fit=exp(p[4])
  E2_fit=exp(p[5])
  alpha1_fit=p[6]
  alpha2_fit=p[7]

  df=n_elements(x)-total(fita)
  chisq=general_tpl_chisq(lnx,lny,x_err=lnx_err,y_err=lny_err,p=p,numerical=numerical)
  reduced_chisq=chisq/df

  C=2*H
  if keyword_set(normalization) then begin
    C=C*reduced_chisq
  endif
  covar=C*(fita#transpose(fita))
  a_std=sqrt(covar[0,0])
  beta1_std=sqrt(covar[1,1])
  beta2_std=sqrt(covar[2,2])
  beta3_std=sqrt(covar[3,3])
  E1_std=sqrt(covar[4,4])*E1_fit
  E2_std=sqrt(covar[5,5])*E2_fit
  alpha1_std=sqrt(covar[6,6])
  alpha2_std=sqrt(covar[7,7])

  print,'**********GT fitting results**********'
  print,'A='+strcompress(string(a_fit),/rem)+'±'+strcompress(string(a_std))
  print,'beta1='+strcompress(string(beta1_fit),/rem)+'±'+strcompress(string(beta1_std))
  print,'beta2='+strcompress(string(beta2_fit),/rem)+'±'+strcompress(string(beta2_std))
  print,'beta3='+strcompress(string(beta3_fit),/rem)+'±'+strcompress(string(beta3_std))
  print,'E1='+strcompress(string(E1_fit),/rem)+'±'+strcompress(string(E1_std))
  print,'E2='+strcompress(string(E2_fit),/rem)+'±'+strcompress(string(E2_std))
  print,'alpha1='+strcompress(string(alpha1_fit),/rem)+'±'+strcompress(string(alpha1_std))
  print,'alpha2='+strcompress(string(alpha2_fit),/rem)+'±'+strcompress(string(alpha2_std))

  ;; print different limit
  case limit of
    2:begin
    print,'---------Double power-law limit:----------'
    print,'A='+strcompress(string(a_fit),/rem)+'±'+strcompress(string(a_std))
    print,'beta1='+strcompress(string(beta2_fit),/rem)+'±'+strcompress(string(beta2_std))
    print,'beta2='+strcompress(string(beta3_fit),/rem)+'±'+strcompress(string(beta3_std))
    print,'E0='+strcompress(string(E2_fit),/rem)+'±'+strcompress(string(E2_std))
    print,'alpha='+strcompress(string(alpha2_fit),/rem)+'±'+strcompress(string(alpha2_std))
  end
  3:begin
  gamma1_fit=beta1_fit+(beta2_fit-beta1_fit)/(1+(E1_fit/min(x))^alpha1_fit)+(beta3_fit-beta2_fit)/(1+(E2_fit/min(x))^alpha2_fit)
  gamma2_fit=beta1_fit+(beta2_fit-beta1_fit)/(1+(E1_fit/max(x))^alpha1_fit)+(beta3_fit-beta2_fit)/(1+(E2_fit/max(x))^alpha2_fit)
  dgamma1_dp=[0,1 - 1d/((E1_fit/min(x))^alpha1_fit + 1),$
    1d/((E1_fit/min(x))^alpha1_fit + 1) - 1d/((E2_fit/min(x))^alpha2_fit + 1),1d/((E2_fit/min(x))^alpha2_fit + 1),$
    0,0,0,0]
  dgamma2_dp=[0,1 - 1d/((exp(p[4])/max(x))^p[6] + 1),$
    1d/((E1_fit/max(x))^alpha1_fit + 1) - 1d/((E2_fit/max(x))^alpha2_fit + 1),1d/((E2_fit/max(x))^alpha2_fit + 1),$
    0,0,0,0]
  gamma1_std=sqrt(total(dgamma1_dp#covar#dgamma1_dp))
  gamma2_std=sqrt(total(dgamma2_dp#covar#dgamma2_dp))
  print,'--------log parabola limit:-------'
  print,'gamma1='+strcompress(string(gamma1_fit),/rem)+'±'+strcompress(string(gamma1_std))
  print,'gamma2='+strcompress(string(gamma2_fit),/rem)+'±'+strcompress(string(gamma2_std))
end
4:begin
Ec_fit=E2_fit*(alpha2_fit/(beta3_fit-beta2_fit))^(1/alpha2_fit)
dlogEc_dp=[0,0,-1d/(p[7]*(p[2] - p[3])),1d/(p[7]*(p[2] - p[3])),0,1,0,1d/p[7]^2d - alog(-p[7]/(p[2] - p[3]))/p[7]^2d]
Ec_std=Ec_fit*sqrt(total(dlogEc_dp#covar#dlogEc_dp))
print,'---------Exponential cut-off limit:----------'
print,'beta='+strcompress(string(beta2_fit),/rem)+'±'+strcompress(string(beta2_std))
print,'Ec='+strcompress(string(Ec_fit),/rem)+'±'+strcompress(string(Ec_std))
print,'alpha='+strcompress(string(alpha2_fit),/rem)+'±'+strcompress(string(alpha2_std))
end
5:begin
w0_fit=E1_fit/beta2_fit
dw0_dp=[0,0,-exp(p[4])/(p[2]^2),0,1d/p[2],0,0,0]
w0_std=sqrt(total(dw0_dp#covar#dw0_dp))
print,'--------Kappa(Generalized Lorentzian) limit:----------'
print,'W0='+strcompress(string(w0_fit),/rem)+'±'+strcompress(string(w0_std),/rem)
end
else:
endcase

if keyword_set(normalization) then begin
  chisq=df
  reduced_chisq=1d
end
if chisq gt 0d and chisq lt 1d+3 and k lt itmax then begin
  Q=1-igamma(0.5*df,0.5*chisq)
  flux_fit=exp(general_tpl(lnx,p=p))

  ECVI=reduced_chisq/(n_elements(x)-1)+2d*total(fita)/double(n_elements(x))
  print,'*********goodness of fit**********'
  print,'chisq='+strcompress(string(chisq),/rem)
  print,'reduced_chisq='+strcompress(string(reduced_chisq),/rem)
  print,'Q='+strcompress(string(Q),/rem)
  print,'ECVI='+strcompress(string(ECVI),/rem)
endif else begin
  print, '-------------failed to converge-----------'
  chisq = 1d+5
endelse
end

function general_tpl,x,p=p
  f=p[0]-p[1]*x+(p[1]-p[2])/p[6]*alog(1+exp(p[6]*(x-p[4])))+(p[2]-p[3])/p[7]*alog(1+exp(p[7]*(x-p[5])))
  return,f
end

function general_tpl_chisq,x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical
  f=general_tpl(x,p=p)
  df_dx=-p[1]+(p[1]-p[2])/(1+exp(-p[6]*(x-p[4])))+(p[2]-p[3])/(1+exp(-p[7]*(x-p[5])))

  if keyword_set(numerical) then begin
    df_dx=deriv(x,y) ;three-point (quadratic) Lagrangian interpolation to compute the derivative
  endif

  return,total((f-y)^2/(y_err^2+(df_dx*x_err)^2))
end

function gt_dJ_dp,x,p=p
  f=general_tpl(x,p=p)

  d1=x-p[4]
  d2=x-p[5]
  df0=replicate(1d,n_elements(x))
  df1=-x+alog(1+exp(p[6]*d1))/p[6]
  df2=-alog(1+exp(p[6]*d1))/p[6]+alog(1+exp(p[7]*d2))/p[7]
  df3=-alog(1+exp(p[7]*d2))/p[7]
  df4=-(p[1]-p[2])/(1+exp(-p[6]*d1))
  df5=-(p[2]-p[3])/(1+exp(-p[7]*d2))
  df6=(p[1]-p[2])/p[6]*(d1/(1+exp(-p[6]*d1))-alog(1+exp(p[6]*d1))/p[6])
  df7=(p[2]-p[3])/p[7]*(d2/(1+exp(-p[7]*d2))-alog(1+exp(p[7]*d2))/p[7])

  return,[[df0],[df1],[df2],[df3],[df4],[df5],[df6],[df7]]
end

function gt_dchisq_dp,x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical
  f=general_tpl(x,p=p)
  chisq=general_tpl_chisq(x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical)
  df_dp=gt_dJ_dp(x,p=p)
  d1=x-p[4]
  d2=x-p[5]
  df_dx=-p[1]+(p[1]-p[2])/(1+exp(-p[6]*(x-p[4])))+(p[2]-p[3])/(1+exp(-p[7]*(x-p[5])))

  d2f_dxdp0=replicate(0d,n_elements(x))
  d2f_dxdp1=-1+1/(1+exp(-p[6]*d1))
  d2f_dxdp2=-1/(1+exp(-p[6]*d1))+1/(1+exp(-p[7]*d2))
  d2f_dxdp3=-1/(1+exp(-p[7]*d2))
  d2f_dxdp4=-(p[1]-p[2])/(1+exp(-p[6]*d1))^2*exp(-p[6]*d1)*p[6]
  d2f_dxdp5=-(p[2]-p[3])/(1+exp(-p[7]*d2))^2*exp(-p[7]*d2)*p[7]
  d2f_dxdp6=(p[1]-p[2])/(1+exp(-p[6]*d1))^2*exp(-p[6]*d1)*d1
  d2f_dxdp7=(p[2]-p[3])/(1+exp(-p[7]*d2))^2*exp(-p[7]*d2)*d2

  if keyword_set(numerical) then begin
    df_dx=deriv(x,y)
    d2f_dxdp0=0d
    d2f_dxdp1=0d
    d2f_dxdp2=0d
    d2f_dxdp3=0d
    d2f_dxdp4=0d
    d2f_dxdp5=0d
    d2f_dxdp6=0d
    d2f_dxdp7=0d
  endif

  dchisq0=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,0]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp0
  dchisq1=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,1]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp1
  dchisq2=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,2]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp2
  dchisq3=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,3]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp3
  dchisq4=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,4]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp4
  dchisq5=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,5]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp5
  dchisq6=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,6]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp6
  dchisq7=2*(f-y)/(y_err^2+(df_dx*x_err)^2)*df_dp[*,7]-(f-y)^2/((y_err^2+(df_dx*x_err)^2))^2*2*x_err^2*df_dx*d2f_dxdp7

  dchisq=total([[dchisq0],[dchisq1],[dchisq2],[dchisq3],[dchisq4],[dchisq5],[dchisq6],[dchisq7]],1)
  return,dchisq
end

function Hessian_matrix_gt_chisq,x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical

  dp=replicate(1d-4,8)
  dchisq=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical)
  H=dblarr(8,8)

  for i=0,7 do begin
    s=dblarr(8)
    s[i]=1d
    dchisq0=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p-dp*s,numerical=numerical)
    dchisq1=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p+dp*s,numerical=numerical)
    H[*,i]=(dchisq1-dchisq0)/2/dp[i]
  endfor
  return,(H+transpose(H))/2

end

function gt_alpha_seek,x,y,x_err=x_err,y_err=y_err,p=p,dp=dp,numerical=numerical
  ;find the alpha where chisq(p+alpha*dp) reaches minimum
  alpha1=0d
  g1=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p+alpha1*dp,numerical=numerical)
  s1=total(g1*dp)

  H=Hessian_matrix_gt_chisq(x,y,x_err=x_err,y_err=y_err,p=p,numerical=numerical)
  alpha2=1d-2*norm(g1)/norm(H#dp)
  g2=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p+alpha2*dp,numerical=numerical)
  s2=total(g2*dp)
  while (alpha2-alpha1)/alpha2 gt 1d-6 do begin
    if s2 lt 0 then begin
      alpha1=alpha2
      alpha2=2.0*alpha2
      s1=s2
      g2=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p+alpha2*dp)
      s2=total(g2*dp)
    endif else begin
      g=gt_dchisq_dp(x,y,x_err=x_err,y_err=y_err,p=p+(alpha1+alpha2)/2*dp)
      s=total(g*dp)
      if s gt 0 then begin
        alpha2=(alpha1+alpha2)/2
        s2=s
      endif else begin
        alpha1=(alpha1+alpha2)/2
        s1=s
      endelse
    endelse
  endwhile
  return,(alpha1+alpha2)/2d
end
