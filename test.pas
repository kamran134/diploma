uses math;

var
i:integer;
mas:array[1..4] of integer;

procedure massplus(var A:array of integer; n:integer);
var
i:integer;
begin
	for i:=1 to n do
	begin
		a[i]:=a[i]+1;
	end;

end;

BEGIN
	for i:=1 to 4 do
	begin
		write('input mas[',i,']: ');
		readln(mas[i]);
	end;

	massplus(mas, 4);

	for i:=1 to 4 do write(mas[i],' ');
END.
