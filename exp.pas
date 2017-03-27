//<a,b> = a1*b1+a2*b2+a3*b3
//<a,b> = |a|*|b|*cos(alpha)
//alpha = arccos(<a,b>/(|a|*|b|))
Uses math;

var
a1, a2, a3, b1, b2, b3: real;
scalar, modul_ab: real;
alpha: real;

begin

write('input a: ');
read(a1, a2, a3);
write('input b: ');
read(b1, b2, b3);

scalar :=a1*b1+a2*b2+a3*b3;
modul_ab := sqrt((a1*a1+a2*a2+a3*a3)*(b1*b1+b2*b2+b3*b3));
alpha := arccos(scalar/modul_ab);

writeln('angle: ', alpha:8:5, ' radians');
writeln('angle: ', (alpha*180/Pi):8:5, ' gradus');

end.
