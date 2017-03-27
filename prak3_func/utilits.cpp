void matrix_print(double **pro, int N) {
	int i, j;
	cout << "\n--------------------\n";
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			cout << pro[i][j] << " ";
		}
		cout << endl;
	}
	cout << "--------------------\n";
}

void array_print_c(double *mas, int N) {
	int i;
	cout << "\n------------------------\n";
	for(i=0; i<N; i++) printf("array[%d] = %.65f\n", i, mas[i]);
	cout << "-------------------------\n";
}

void array_print_l(double *mas, int N) {
	int i;
	cout << "\n------------------------\n";
	for(i=0; i<N; i++) printf("%.45f\t", mas[i]);
	cout << "-------------------------\n";
}

