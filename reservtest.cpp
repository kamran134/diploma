		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=1; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[0] = ddopri5(6,fcn,0,y,beta[3],eps[2],1.0e0,0.5e0);
		X[0] = y[0];
		Y[0] = y[1];
		V[0] = y[2];
		PSI_V[0] = y[5];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[2], GLOBAL_ERROR[0]);
		cout << "eps = " << eps[2] << " : " << endl << "\t x3: " << X[0]  << endl << "\t y3: " << Y[0] << endl << "\t v3: " << V[0] << endl << "\t psi_v_3: " << PSI_V[0] << endl;
		
