Uses math;

var
x,y:real;
res1, res2:real;

Function f(a,b,c:real):real;
begin
	f:=(2*a-b-sin(c))/(5+abs(c));
end;

begin
write('input x, y: ');
read(x, y);

res1:=f(x,-2*y,1.17)+f(2.2,x,x-y);
res2:=tan(f(x+y,x*y,y-x))+f(3.1,1.4,y-sin(x));

writeln('1: ', res1:8:5);
writeln('2: ', res2:8:5);
end.
