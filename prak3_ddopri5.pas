Uses math;

type
  matrix3 = array[1..4, 1..4] of real;

var
y:array[1..6] of real;
beta:array[1..4] of real;
//i:integer;
k1,k2,k3,k4,k5,y1:array[1..11] of real;
b_cpy,x_nq,b_new,r:array[1..4] of real;
eps2:real=2.e-15;
pro:matrix3;
dbeta:array[1..4] of real;


Function norm(vector:array of real; n:integer):real;
var
	i:integer;
	cur,max_w,sum2:real;
begin
	
	max_w:=0.0;
	
	for i:=0 to n do
	begin
		cur:=abs(vector[i]);
		if cur>max_w then max_w:=cur;
	end;
	
	if max_w=0.0 then norm:=0;
	sum2:=0.0;
	for i:=0 to n do
	begin
		cur:=vector[i]/max_w;
		sum2:=sum2+cur*cur;
	end;
	norm:=max_w*sqrt(sum2);
end;

Procedure fcn(var x:real; y,f:array of real);
var
f_right:array[1..6] of real;
RHO,psi_xv,psi_yvpsi_v:real;
begin
	psi_xv:=y[4]*y[4]*y[3]*y[3];
    psi_yvpsi_v:=(y[5]*y[3]-y[6])*(y[5]*y[3]-y[6]);
    RHO:=sqrt(psi_xv+psi_yvpsi_v);
    
    //правая часть
    //f(v) = v^n, 	n = 1, 	=>	 f'(v) = (v)' = 1;
    f_right[1]:=y[3]*y[3]*y[4]/RHO;									//x^{\cdot}
	f_right[2]:=y[3]*(y[5]*y[3]-y[6])/RHO;								//y^{\cdot}
	f_right[3]:= -y[3]-(y[5]*y[3]-y[6])/RHO; 							//v^{\cdot}
	f_right[4]:=0;														//psi_x^{\cdot}
	f_right[5]:=0;														//psi_y^{\cdot}
	f_right[6]:=y[6]-(y[3]*y[4]*y[4]-y[3]*y[5]*y[5]-y[5]*y[6])/RHO;	//psi_v^{\cdot}
end;

//МЕТОД РУНГЕ-КУТТЫ

Function dsign(a,b:real):real;
begin 
	if b<0 then 
	begin
		a:=abs(a)*(-1.0); 
		dsign:=a;
	end;
	if b>0 then dsign:=abs(a);
	dsign:=0.0;
end;

function max_root(x, lambda:real):real;
begin
	max_root:=max(-0.5+0.5*lambda*x, 0.5-0.5*lambda*x);
end;

	
function ddopri5(n:integer; x:real; y:array of real; xend,eps,hmax,h:real):real;
var
    xmult:real;
    reject:boolean;
    xph,err,denom,fac,hnew,posneg,omega:real;
    nmax:integer=30000;
    i:integer;
    uround:real=2.2205e-16;
    gerror:real=0;
    naccpt:integer=0;
    nrejct:integer=0;
    nfcn:integer=1;
    nstep:integer=0;
begin
	writeln('DDOPRI5');
    posneg:=dsign(1.e0,xend-x);
    hmax:=abs(hmax);
    h:=min(max(1.e-4,abs(h)),hmax);
    h:=dsign(h,posneg);
    reject:=false;
    
    fcn(x,y,k1);
    while 1=1 do
    begin
		if nstep>nmax then break;
		if (x-xend)*posneg+uround>0.e0 then break;
		if (x+h-xend)*posneg > 0.e0 then h:=xend-x;
		nstep:=nstep+1;
		for i:=1 to n do y1[i]:=y[i]+h*0.2e0*k1[i];
		xmult:=x+h*0.2;
		fcn(xmult,y1,k2);
		for i:=1 to n do y1[i]:=y[i]+h*((3.e0/40.e0)*k1[i]+(9.e0/40.e0)*k2[i]);
		xmult:=x+h*0.3e0;
		fcn(xmult,y1,k3);
		for i:=1 to n do y1[i]:=y[i]+h*((44.e0/45.e0)*k1[i]-(56.e0/15.e0)*k2[i]+(32.e0/9.e0)*k3[i]);
		xmult:=x+h*0.8e0;
		fcn(xmult,y1,k4);
		for i:=1 to n do y1[i]:=y[i]+h*((19372.e0/6561.e0)*k1[i]-(25360.e0/2187.e0)*k2[i]+(64448.e0/6561.e0)*k3[i]-(212.e0/729.e0)*k4[i]);
		xmult:=x+h*(8.e0/9.e0);
		fcn(xmult,y1,k5);
		for i:=1 to n do y1[i]:=y[i]+h*((9017.e0/3168.e0)*k1[i]-(355.e0/33.e0)*k2[i]+(46732.e0/5247.e0)*k3[i]+(49.e0/176.e0)*k4[i]-(5103.e0/18656.e0)*k5[i]);
		xph:=x+h;
		fcn(xph,y1,k2);
		for i:=1 to n do y1[i]:=y[i]+h*((35.e0/384.e0)*k1[i]+(500.e0/1113.e0)*k3[i]+(125.e0/192.e0)*k4[i]-(2187.e0/6784.e0)*k5[i]+(11.e0/84.e0)*k2[i]);
		for i:=1 to n do k2[i]:=(71.e0/57600.e0)*k1[i]-(71.e0/16695.e0)*k3[i]+(71.e0/1920.e0)*k4[i]-(17253.e0/339200.e0)*k5[i]+(22.e0/525.e0)*k2[i];
		fcn(xph,y1,k3);
		for i:=1 to n do k4[i]:=(k2[i]-(1.e0/40.e0)*k3[i])*h;
		nfcn:=nfcn+6;
		err:=0;
		for i:=1 to n do
		begin
			denom:=max(1.e-5,max(abs(y1[i]),max(abs(y[i]),2.e0*uround/eps)));
			err:=err+power(k4[i]/denom,2);
		end;
		err:=sqrt(err/double(n));
		fac:=max(0.1e0, min( 5.e0, power(err/eps,0.2e0 )/0.9e0) );
		hnew:=h/fac;
		if err<=eps then
		begin
			naccpt:=naccpt+1;
			for i:=1 to n do
			begin
				k1[i]:=k3[i];
				y[i]:=y1[i];
			end;
			x:=xph;
			if abs(hnew)>hmax then hnew:=posneg*hmax;
			if reject=true then 
			begin
				hnew:=posneg*min(abs(hnew),abs(h)); 
				reject:=false;
			end
			else reject:=true;
			if naccpt >= 1 then nrejct:=nrejct+1;
			omega:=max_root(x, y[3]);
			gerror:=err + gerror*power(exp(1), h*omega); //gerror это то, что в теории называется лямбдой маленькой, а err это лямбда большая (треугольная)
		end;
		h:=hnew;
		
		
		writeln(x);
		writeln('     x(',x:8:16,') = ', y[0]);
		writeln('     y(',x:8:16,') = ', y[1]);
		writeln('     v(',x:8:16,') = ', y[2]);
		writeln('     psi_x(',x:8:16,') = ', y[3]);
		writeln('     psi_y(',x:8:16,') = ', y[4]);
		writeln('     psi_v(',x:8:16,') = ', y[5]);
		writeln('-------------------------------');
		
    end;
    ddopri5:=gerror;
end;



Procedure f(var beta,res:array of real);
var
//y:array[1..6] of real;
psi_xv,psi_yvpsi_v,RHO:real;
begin
	y[1]:=0; //x_0
	y[2]:=10; //y_0
	y[3]:=0; //v_0
	y[4]:=beta[1];
	y[5]:=beta[2];              
	y[6]:=beta[3];
	ddopri5(6,0,y,beta[4],1.e-11,1.0e0,0.5e0);
	psi_xv:=y[4]*y[4]*y[3]*y[3];
    psi_yvpsi_v:=(y[5]*y[3]-y[6])*(y[5]*y[3]-y[6]);
    RHO:=sqrt(psi_xv+psi_yvpsi_v); //RHO(T)
    
	res[1]:=y[1]-0.1386215299798196;//+0.0000000000000000277555756156289135105907917022705078125;//-0.1386315299798196; //x_T
	res[2]:=y[2]-8.8869461060126618;//-8.8869461060126618; //y_T
	res[3]:=y[3]-0.8583067592798990;//-0.8583067592798990; //v_T
	res[4]:=RHO-res[3]*y[3]-1;
end;

Procedure gradf(beta:array of real; pro:matrix3);
var
i,j:integer;
N:integer=4;
beta_p, beta_m, res_p, res_m:array[0..3] of real;
h:real=1.e-6;
begin
	for j:=1 to n do
	begin 
		for i:=1 to n do
		begin
			beta_p[i]:=beta[i];
			beta_m[i]:=beta[i];
		end;
		
		beta_p[j]:=beta_p[j]+h;
		beta_m[j]:=beta_m[j]-h;
		
		f(beta_p,res_p);
		f(beta_m,res_m);
		
		for i:=1 to n do pro[i][j]:=(res_p[i]-res_m[i])/(2*h);
	end;
end;

//gauss
procedure swap_el(a,b:real);
var tmp:real;
begin
tmp:=a;
a:=b;
b:=tmp;
end;

procedure masswap(a:matrix3; b:array of real; n,i,j:integer);
var
k:integer;
begin
    for k:=0 to n-1 do swap_el(a[i][k], a[j][k]);
    swap_el(b[i], b[j]);
end;

procedure gauss(var a:matrix3;x,b:array of real; n:integer);
var
i, j, k, mxi:integer;
mx:real; //max
begin
    for i:=1 to n do
    begin
        mx:=a[i][i];
        mxi:=i;
        for j:=i+1 to n-1 do
            if a[j][i]>mx then
            begin
                mx:=a[j][i];
                mxi:=j;
            end;
        masswap(a,b,n,i,mxi);
        for k:=i+1 to n-1 do
        begin
            for j:=i+1 to n-1 do
            begin
                if (a[k][i]>0.0) or (a[k][i]<0.0) then a[k][j]:=a[k][j]/a[k][i]-a[i][j]/a[i][i];
            end;
            if (a[k][i]>0.0) or (a[k][i]<0.0) then b[k]:=b[k]/a[k][i]-b[i]/a[i][i];
            a[k][i]:=0.0;
        end;
        b[i]:=b[i]/a[i][i];
        for j:=i+1 to n-1 do a[i][j]:=a[i][j]/a[i][i];
        a[i][i]:=1.0;
    end;
    for i:=n-1 downto 0 do
    begin
        for j:=n-1 downto i+1 do
        begin
            b[i]:=b[i]-a[i][j]*x[j];
        end;
        x[i]:=b[i];
    end;
end;

procedure matrix_multiplication(var a:matrix3; x,b:array of real;n:integer);
var i,j:integer;
begin
    for i:=1 to n do
    begin
        b[i]:= 0;
        for j:=1 to n do
        begin
            b[i]:=b[i]+(a[i][j]*x[j]);
        end;
    end;
end;

procedure linear_sys(var a:matrix3; x,b:array of real; n:integer);
var
a_cpy:matrix3;

i,j:integer;
begin
    for i:= 0 to n-1 do
    begin 
		for j:=1 to n do a_cpy[i][j]:= a[i][j];
	end; //копирование матрицы
    
    for i:=1 to n do b_cpy[i]:= b[i]; //копирование вектора
    
    gauss(a_cpy, x_nq, b_cpy, n);
    matrix_multiplication(a, x_nq, b_new, n);
    
    for i:=1 to n do 
    begin
		for j:=1 to n do a_cpy[i][j]:=a[i][j];
	end;
	
    for i:=1 to n do b_cpy[i]:=b[i]-b_new[i];
    
    gauss(a_cpy, r, b_cpy, n);
    for i:= 0 to n-1 do x[i]:=x_nq[i]-r[i];
end;





function NEWTON(beta:array of real):integer;
	var
	i, k, flag:integer;
	N:integer=4;
	res, res_w, beta_w, res_m:array[0..3] of real;
	gamma:real;
begin
	f(beta, res);
	
	for k:=0 to 14 do
	begin
		writeln( '-------');
		writeln('k = ',k);
		writeln( '-------');
		if (abs(res[0])<eps2) and (abs(res[1])<eps2) and (abs(res[2])<eps2) then 
		begin
			writeln('Ended by ',k,' iteration');
			NEWTON:=k;
		end;
		gradf(beta, pro);
		//kramer(pro, dbeta, res);
		for i:=1 to n do res_m[i]:=-res[i];
		linear_sys(pro, dbeta, res_m, 4);
		
		
		gamma:=1.0;
		flag:=-1;
		
		while gamma>sqrt(eps2) do
		begin
			for i:=1 to n do beta_w[i]:=beta[i]+gamma*dbeta[i];
			f(beta_w, res_w);
			if norm(res_w,N)<norm(res,N) then
			begin
				flag:=0; 
				break;
			end;
			gamma:=gamma/2.0;
		end;
		if flag=-1 then
		begin
			writeln('Broken on ',k,'`s iteration');
			NEWTON:=-1;
		end;
		for i:=1 to n do
		begin
			beta[i]:=beta_w[i];
			res[i]:=res_w[i];
		end;
	end;
writeln('not enough iteration');
NEWTON:=-2;
end;

begin
	beta[1]:=1;
	beta[2]:=1;
	beta[3]:=1;
	beta[4]:=2;
	
	NEWTON(beta);
end.
