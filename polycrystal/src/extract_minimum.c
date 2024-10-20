typedef struct{
	double value;
	int index;
}MinResult;

MinResult minimumArray(double* array, int size){
	MinResult result;
	result.value = array[0];
	result.index = 0;

	for (int i = 1; i < size; i++){
		if (array[i] < result.value) {
			result.value = array[i];
			result.index = i;
		}
	}
	return result;
}