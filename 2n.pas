var
n,m:integer;
two:integer;
begin
write('input m: ');
readln(m);

n:=-1;
two:=1;

while two<m do
	begin
		two:=two*2;
		n:=n+1;
	end;
if m<=0 then
	begin
		writeln('n not exists');
	end
else
	begin
		writeln('max{n} = ', n);
		writeln('2^', n, ' < ', m);
	end;
end.
