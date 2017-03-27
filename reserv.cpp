//ddopri5(6,fcn,0,y,1.0e0,eps[2],1.0e0,0.5e0);

//В точке t = 0.0260000000000000
	cout << "t = 0.26" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=1; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,0.026,eps[i],1.0e0,0.5e0);
		X[i] = y[0];
		Y[i] = y[1];
		V[i] = y[2];
		PSI_V[i] = y[5];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[0] << " : " << endl << "\t x1: " << X[0]  << endl << "\t y1: " << Y[0] << endl << "\t v1: " << V[0] << endl << "\t psi_v_1: " << PSI_V[0] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[1] << " : " << endl << "\t x2: " << X[1]  << endl << "\t y2: " << Y[1] << endl << "\t v2: " << V[1] << endl << "\t psi_v_2: " << PSI_V[1] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[2] << " : " << endl << "\t x3: " << X[2]  << endl << "\t y3: " << Y[2] << endl << "\t v3: " << V[2] << endl << "\t psi_v_3: " << PSI_V[2] << endl;
	cout << "**********************" << endl;
	
	x1 = X[0] - X[1];
	x2 = X[1] - X[2];
	x3 = x1/x2;
	y1 = Y[0] - Y[1];
	y2 = Y[1] - Y[2];
	y3 = y1/y2;
	v1 = V[0] - V[1];
	v2 = V[1] - V[2];
	v3 = v1/v2;
	psi_v1 = PSI_V[0] - PSI_V[1];
	psi_v2 = PSI_V[1] - PSI_V[2];
	psi_v3 = psi_v1/psi_v2;
	
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("x1 - x2 = %.16f\nx2 - x3 = %.16f\n(x1-x2)/(x2-x3) = %.16f\n",x1,x2,x3);
	printf("\ny1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",y1,y2,y3);
	printf("\nv1 - v2 = %.16f\nv2 - v3 = %.16f\n(v1-v2)/(v2-v3) = %.16f\n",v1,v2,v3);
	printf("\npsi_v1 - psi_v2 = %.16f\npsi_v2 - psi_v3 = %.16f\n(psi_v1-psi_v2)/(psi_v2-psi_v3) = %.16f\n",psi_v1,psi_v2,psi_v3);
	printf("\nd1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
	
	//В точке t = 1/2 T
	cout << "t = 1/2 T" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=1; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,beta[3]/2.,eps[i],1.0e0,0.5e0);
		X[i] = y[0];
		Y[i] = y[1];
		V[i] = y[2];
		PSI_V[i] = y[5];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[0] << " : " << endl << "\t x1: " << X[0]  << endl << "\t y1: " << Y[0] << endl << "\t v1: " << V[0] << endl << "\t psi_v_1: " << PSI_V[0] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[1] << " : " << endl << "\t x2: " << X[1]  << endl << "\t y2: " << Y[1] << endl << "\t v2: " << V[1] << endl << "\t psi_v_2: " << PSI_V[1] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[2] << " : " << endl << "\t x3: " << X[2]  << endl << "\t y3: " << Y[2] << endl << "\t v3: " << V[2] << endl << "\t psi_v_3: " << PSI_V[2] << endl;;
	cout << "**********************" << endl;
	
	x1 = X[0] - X[1];
	x2 = X[1] - X[2];
	x3 = x1/x2;
	y1 = Y[0] - Y[1];
	y2 = Y[1] - Y[2];
	y3 = y1/y2;
	v1 = V[0] - V[1];
	v2 = V[1] - V[2];
	v3 = v1/v2;
	psi_v1 = PSI_V[0] - PSI_V[1];
	psi_v2 = PSI_V[1] - PSI_V[2];
	psi_v3 = psi_v1/psi_v2;
	
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("x1 - x2 = %.16f\nx2 - x3 = %.16f\n(x1-x2)/(x2-x3) = %.16f\n",x1,x2,x3);
	printf("\ny1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",y1,y2,y3);
	printf("\nv1 - v2 = %.16f\nv2 - v3 = %.16f\n(v1-v2)/(v2-v3) = %.16f\n",v1,v2,v3);
	printf("\npsi_v1 - psi_v2 = %.16f\npsi_v2 - psi_v3 = %.16f\n(psi_v1-psi_v2)/(psi_v2-psi_v3) = %.16f\n",psi_v1,psi_v2,psi_v3);
	printf("\nd1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
	
	//В точке t = T
	cout << "t = T" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=1; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,beta[3],eps[i],1.0e0,0.5e0);
		X[i] = y[0];
		Y[i] = y[1];
		V[i] = y[2];
		PSI_V[i] = y[5];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[0] << " : " << endl << "\t x1: " << X[0]  << endl << "\t y1: " << Y[0] << endl << "\t v1: " << V[0] << endl << "\t psi_v_1: " << PSI_V[0] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[1] << " : " << endl << "\t x2: " << X[1]  << endl << "\t y2: " << Y[1] << endl << "\t v2: " << V[1] << endl << "\t psi_v_2: " << PSI_V[1] << endl;
	cout <<fixed<<setprecision(16)<< "eps = " << eps[2] << " : " << endl << "\t x3: " << X[2]  << endl << "\t y3: " << Y[2] << endl << "\t v3: " << V[2] << endl << "\t psi_v_3: " << PSI_V[2] << endl;;
	cout << "**********************" << endl;
	
	x1 = X[0] - X[1];
	x2 = X[1] - X[2];
	x3 = x1/x2;
	y1 = Y[0] - Y[1];
	y2 = Y[1] - Y[2];
	y3 = y1/y2;
	v1 = V[0] - V[1];
	v2 = V[1] - V[2];
	v3 = v1/v2;
	psi_v1 = PSI_V[0] - PSI_V[1];
	psi_v2 = PSI_V[1] - PSI_V[2];
	psi_v3 = psi_v1/psi_v2;
	
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("x1 - x2 = %.16f\nx2 - x3 = %.16f\n(x1-x2)/(x2-x3) = %.16f\n",x1,x2,x3);
	printf("\ny1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",y1,y2,y3);
	printf("\nv1 - v2 = %.16f\nv2 - v3 = %.16f\n(v1-v2)/(v2-v3) = %.16f\n",v1,v2,v3);
	printf("\npsi_v1 - psi_v2 = %.16f\npsi_v2 - psi_v3 = %.16f\n(psi_v1-psi_v2)/(psi_v2-psi_v3) = %.16f\n",psi_v1,psi_v2,psi_v3);
	printf("\nd1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
