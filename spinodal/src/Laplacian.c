void Laplacian(double* array, double* laplacian){
	
	double up,down,left,right,center;

	for(int i = 0; i < N; i++){

		if(i >= H){up = array[i-H];}
		else{up = array[N-H+i];}

		if(N-H > i){down = array[i+H];}
		else{down = array[i+H-N];}

		if(i%W == 0){left = array[i+W-1];}
		else{left = array[i-1];}

		if((i+1)%W == 0){right = array[i+1-W];}
		else{right = array[i+1];}

		center = array[i];

		laplacian[i] = (right+left-2*center)/(step*step) + (up+down-2*center)/(step*step);}

}