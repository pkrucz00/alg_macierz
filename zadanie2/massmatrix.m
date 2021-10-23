function [A] = massmatrix(riga,nxx,pxx,rxx)

%riga: 0=iga, 1=mes, 2=riga

%funkcja liczaca ilosc funkcji bazowych
compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1;

if riga==0 
  %uniform knot vector B-splines C^(p-1)
  knot2 = simple_knot(nxx,pxx);
end
if riga==1
  %Lagrange
  knot2(1:pxx+1)=0; 
  k=pxx+2;
  for i=1:nxx
    knot2(k:k+1)=1/nxx*i;
    k=k+2;
  end 
  knot2(2*nxx+pxx+1:2*nxx+2*(pxx+1)-2)=1;
end
if riga==2
  if rxx==0
    warning('For riga parameter set to "2" the rxx parameter cannot be set to "0". Please check parameters passed to "masslowrank" function. The function will behave as if the riga parameter was set to "0"')
  end
  
  knot2(1:pxx+1)=0;
  k=pxx+2;
  
  for i=1:nxx
    knot2(k)=1/nxx*i;
    k=k+1;
    
    if mod(i,rxx)==0 && i~=nxx
      knot2(k)=1/nxx*i;
      k=k+1;
    end
  end

  knot2(k:k+pxx-1)=1;
end
knot2

knot_vectorx=knot2;
knot_vectory=knot2;  

px = compute_p(knot_vectorx);
py = compute_p(knot_vectory);

elementsx = number_of_elements(knot_vectorx,px);
elementsy = number_of_elements(knot_vectory,py);

nx = number_of_dofs(knot_vectorx,px);
ny = number_of_dofs(knot_vectory,py);

A = sparse(nx*ny,nx*ny);

%calki z B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
%(i,k=1,...,Nx; j,l=1,...,Ny)
%petla po elementach w osi x
for ex = 1:elementsx;
  % zakres funkcji niezerowych nad elementem
  [xl,xh] = dofs_on_element(knot_vectorx,px,ex);
  % zakres elementu (lewy i prawy brzeg po x)
  [ex_bound_l,ex_bound_h] = element_boundary(knot_vectorx,px,ex);
  %petla po elementach w osi y
  for ey = 1:elementsy
    % zakres funkcji niezerowych nad elementem
    [yl,yh] = dofs_on_element(knot_vectory,py,ey);
    % zakres elementu (lewy i prawy brzeg po x)
    [ey_bound_l,ey_bound_h] = element_boundary(knot_vectory,py,ey);
    % jakobian - rozmiar elementu
    Jx = ex_bound_h - ex_bound_l;
    Jy = ey_bound_h - ey_bound_l;
    J = Jx*Jy;
    % petla po funkcjach niezerowych nad danym elementem
    for bi = xl:xh
      for bj = yl:yh
        for bk = xl:xh
          for bl = yl:yh
          
          % punkty kwadratury w osi x nad elementem
          qpx = quad_points(ex_bound_l,ex_bound_h,2*px+2*py+1);
          % punkty kwadratury w osi y nad elementem
          qpy = quad_points(ey_bound_l,ey_bound_h,2*px+2*py+1);
          % wagi kwadratury w osi x nad elementem
          qwx = quad_weights(ex_bound_l,ex_bound_h,2*px+2*py+1);
          % wagi kwadratury w osi y nad elementem
          qwy = quad_weights(ey_bound_l,ey_bound_h,2*px+2*py+1);
          % petla po punktach kwadratury
          for iqx = 1:size(qpx,2)
            for iqy = 1:size(qpy,2)
              % definicja funkcji ksztaltu
              % B^x_k(x)
              funk= compute_spline(knot_vectorx,px,bi,qpx(iqx));
              % B^y_l(y)
              funl= compute_spline(knot_vectory,py,bj,qpy(iqy));
              % B^x_i(x)
              funi= compute_spline(knot_vectorx,px,bk,qpx(iqx));
              % B^y_j(y)
              funj= compute_spline(knot_vectory,py,bl,qpy(iqy));
              % B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
              fun = funi*funj*funk*funl;
              
              % Calki z B^x_i(x) B^y_j(y) B^x_k(x) B^y_l(y)
              % (i,k=1,...,Nx; j,l=1,...,Ny)
              int = fun*qwx(iqx)*qwy(iqy)*J;
              if (int~=0)
                A((bj-1)*nx+bi,(bl-1)*nx+bk) = A((bj-1)*nx+bi,(bl-1)*nx+bk) + int;
              end
            end
          end
          end
        end
      end
    end
  end
end

spy(A);

return

A11 = A(1:nx*ny/2,1:nx*ny/2);
A12 = A(1:nx*ny/2,nx*ny/2:nx*ny);

%A11
[U,S,V]=svd(A11);

size=nx*ny/2;
for i=1:size
  diag(i)=S(i,i);
end

%diag

figure(1); 
arg=[1:1:size];
set(gca,'Fontsize',60); 
plot(arg,diag,'LineWidth',3);
title('A11');
hold on

%A11inv
A11inv=inverse(A11);
[U,S,V]=svd(A11inv);

size=nx*ny/2;
for i=1:size
  diag(i)=S(i,i);
end

%diag

figure(2); 
arg=[1:1:size];
set(gca,'Fontsize',60); 
plot(arg,diag,'LineWidth',3);
title('A11^{-1}');
hold on

%A12

[U,S,V]=svd(A12);
size=nx*ny/2;
for i=1:size
  diag(i)=S(i,i);
end

%diag

figure(3); 
arg=[1:1:size];
set(gca,'Fontsize',60); 
plot(arg,diag,'LineWidth',3);
title('A12');

%A11^{-11}A12
A11invA12=A11inv*A12;

[U,S,V]=svd(A11invA12);
size=nx*ny/2;
for i=1:size
  diag(i)=S(i,i);
end

%diag

figure(4); 
arg=[1:1:size];
set(gca,'Fontsize',60); 
plot(arg,diag,'LineWidth',3);
title('A11^{-1}*A12');

%funkcja wyliczajaca stopien wielomianow
function p=compute_p(knot_vector)

%pierwszy wpis w knot_vector
  initial = knot_vector(1);
%dlugosc knot_vector
  kvsize = size(knot_vector,2);
  p = 0;

%sprawdzenie ilosci powtorzen pierwszego wpisu w knot_vector  
  while (p+2 <= kvsize) && (initial == knot_vector(p+2))
    p = p+1;
  endwhile
  
  return
end
  
  
  
%funkcja sprawdzajaca poprawnosc knot_vector
function t=check_sanity(knot_vector,p)

  initial = knot_vector(1);
  kvsize = size(knot_vector,2);

  t = true;
  counter = 1;

%jesli ilosc powtorzen na poczatku nie jest zgodna ze stopniem wielomianow
  for i=1:p+1
    if (initial ~= knot_vector(i))
%zwroc falsz
      t = false;
      return
    end
  end

%jesli w srodku knot_vector jest za duzo powtorzen
  for i=p+2:kvsize-p-1
    if (initial == knot_vector(i))
      counter = counter + 1;
      if (counter > p)
%zwroc falsz
        t = false;
        return
      end
    else
      initial = knot_vector(i);
      counter = 1;
    end
  end

  initial = knot_vector(kvsize);
  
%jesli ilosc powtorzen na poczatku nie jest zgodna ze stopniem wielomianow
  for i=kvsize-p:kvsize
    if (initial ~= knot_vector(i))
%zwroc falsz
      t = false;
      return
    end
  end
  
%jesli kolejny element knot_vector jest mniejszy niz poprzedni
  for i=1:kvsize-1
    if (knot_vector(i)>knot_vector(i+1))
%zwroc falsz
      t = false;
      return
    end
  end

  return

end



%funkcja wyliczajaca rekurencyjnie funkcje bazowa zgodnie z wzorem Cox-de-Boor
function y=compute_spline(knot_vector,p,nr,x)
  
%funkcja (x-x_i)/(x_{i-p}-x_i)
  fC= @(x,a,b) (x-a)/(b-a);
%funkcja (x_{i+p+1}-x)/(x_{i+p+1}-x_{i+1})
  fD= @(x,c,d) (d-x)/(d-c);
  
%x_i  
  a = knot_vector(nr);
%x_{i-p} 
  b = knot_vector(nr+p);
%x_{i+1}
  c = knot_vector(nr+1);
%x_{i+p+1}
  d = knot_vector(nr+p+1);

%funkcja liniowa dla p=0
  if (p==0)
    y = 0 .* (x < a) + 1 .* (a <= x & x <= d) + 0 .* (x > d);
    return
  end

%B_{i,p-1}  
  lp = compute_spline(knot_vector,p-1,nr,x);
%B_{i+1,p-1}
  rp = compute_spline(knot_vector,p-1,nr+1,x);
  
%(x-x_i)/(x_{i-p)-x_i)*B_{i,p-1}  
  if (a==b)
%jesli wezly w knot_vector sie powtarzaja trzeba to uwzglednic
    y1 = 0 .* (x < a) + 1 .* (a <= x & x <= b) + 0 .* (x > b);
  else
    y1 = 0 .* (x < a) + fC(x,a,b) .* (a <= x & x <= b) + 0 .* (x > b);
  end
  
%(x_{i+p+1}-x)/(x_{i+p+1)-x_{i+1})*B_{i+1,p-1}
  if (c==d)
%jesli wezly w knot_vector sie powtarzaja trzeba to uwzglednic
    y2 = 0 .* (x < c) + 1 .* (c < x & x <= d) + 0 .* (d < x);
  else
    y2 = 0 .* (x < c) + fD(x,c,d) .* (c < x & x <= d) + 0 .* (d < x);
  end
  
  y = lp .* y1 + rp .* y2;
  return
  
end

  
  
% Ilosc elementow w knot vector
function n = number_of_elements(knot_vector,p)
  
  initial = knot_vector(1);
  kvsize = size(knot_vector,2);
  n = 0;
  
  for i=1:kvsize-1
    if (knot_vector(i) ~= initial)
      initial = knot_vector(i);
      n = n+1;
    end
  end
  
end
  
% Tworzy prosty knot vector bez powotrzen w srodku
function knot = simple_knot(elems, p)
  pad = ones(1, p);
  knot = [0 * pad, 0:elems, elems * pad];
  knot = knot/elems;
end

% Ilosc funkcji bazowych w knot wektorze
function n = number_of_dofs(knot,p)
  n = length(knot) - p - 1;
end


function first = first_dof_on_element(knot_vector,p,elem_number)
 [l,h] = element_boundary(knot_vector,p,elem_number);
 first = lookup(knot_vector, l) - p;
end
  
function [low,high] = element_boundary(knot_vector,p,elem_number)
  initial = knot_vector(1);
  kvsize = size(knot_vector,2);
  k = 0;
  low=0;
  high=0;
  
  for i=1:kvsize
    if (knot_vector(i) ~= initial)
      initial = knot_vector(i);
      k = k+1;
    end
    if (k == elem_number)
      low = knot_vector(i-1);
      high = knot_vector(i);
      return;
    end
  end
end

% Zwraca zakres (indeksy) funckji bedacych niezerowymi na zadanym wektorze wezlow
function [low,high] = dofs_on_element(knot_vector,p,elem_number)
  low = first_dof_on_element(knot_vector,p,elem_number);
  %poniewaz mamy dokladnie p+1 niezerowych funkcji nad elementem
  high = low + p;
end

% Row vector of points of the k-point Gaussian quadrature on [a, b]
function xs = quad_points(a, b, k)
  % mapowanie punktow
  map = @(x) 0.5 * (a * (1 - x) + b * (x + 1));
  switch (k)
    case 1
      xs = [0];
    case 2
      xs = [-1/sqrt(3), ...
             1/sqrt(3)];
    case 3
      xs = [-sqrt(3/5), ...
             0,         ...
             sqrt(3/5)];
    case 4
      xs = [-sqrt((3+2*sqrt(6/5))/7), ...
             sqrt((3-2*sqrt(6/5))/7), ...
             sqrt((3-2*sqrt(6/5))/7), ...
             sqrt((3+2*sqrt(6/5))/7)];
    case 5
      xs = [-1/3*sqrt(5+2*sqrt(10/7)), ...
            -1/3*sqrt(5-2*sqrt(10/7)), ...
             0,                        ...
             1/3*sqrt(5-2*sqrt(10/7)), ...
             1/3*sqrt(5+2*sqrt(10/7))];
    otherwise
      xs = [-1/3*sqrt(5+2*sqrt(10/7)), ...
            -1/3*sqrt(5-2*sqrt(10/7)), ...
             0,                        ...
             1/3*sqrt(5-2*sqrt(10/7)), ...
             1/3*sqrt(5+2*sqrt(10/7))];
  endswitch
  xs = map(xs);
end

% Row vector of weights of the k-point Gaussian quadrature on [a, b]
function ws = quad_weights(a, b, k)
  switch (k)
    case 1
      ws = [2];
    case 2
      ws = [1, 1];
    case 3
      ws = [5/9, ...
            8/9, ...
            5/9];
    case 4
      ws = [(18-sqrt(30))/36, ...
            (18+sqrt(30))/36, ...
            (18+sqrt(30))/36, ...
            (18-sqrt(30))/36];
    case 5
      ws = [(322-13.0*sqrt(70))/900, ...
            (322+13.0*sqrt(70))/900, ...
            128/225,                 ...
            (322+13.0*sqrt(70))/900, ...
            (322-13.0*sqrt(70))/900];
    otherwise
      ws = [(322-13.0*sqrt(70))/900, ...
            (322+13.0*sqrt(70))/900, ...
            128/225,                 ...
            (322+13.0*sqrt(70))/900, ...
            (322-13.0*sqrt(70))/900];
  end
  % Gaussian quadrature is defined on [-1, 1], we use [a, b]
  %ws = ws * (b-a)/2;
end
end