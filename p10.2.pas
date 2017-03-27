var
A:array[1..10] of integer;
sred:real;
i,minus:integer;

begin
	
	for i:=1 to 10 do
	begin
		write('A[', i, '] = ');
		readln(A[i]);
	end;
	
	minus:=0;
	sred:=0;
	for i:=1 to 10 do
	begin
		if A[i]<0 then
		begin
			minus:=minus+1;
		end;
		sred:=sred+A[i];
	end;
	sred:=sred/10;
	
	writeln('menfi elementlerin sayi: ', minus);
	writeln('ededi orta: ', sred:8:16);
end.
