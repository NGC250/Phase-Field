void ImagefromArray(double* array, int rows, int columns, int index, int max_iteration, const char* filename) {
    
	char files[128];
    char buffer[128];

    snprintf(buffer, sizeof(buffer), "%s", filename);
    int num_digits = snprintf(NULL, 0, "%d", max_iteration);

    if (snprintf(files, sizeof(files), "%s%0*d.pgm", buffer, num_digits, index) >= sizeof(files)) {
        fprintf(stderr, "Error: File path buffer too small\n");
        exit(EXIT_FAILURE);
    }
	
    FILE *pgmimg = fopen(files, "wb");
	
    if (pgmimg == NULL){printf("Error: Could not open file %s for writing\n", files);return;}

    fprintf(pgmimg, "P2\n%d %d\n255\n", W, H);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            int temp = (int)(array[i * columns + j] * 255);
            fprintf(pgmimg, "%d ", temp);
        }
        fprintf(pgmimg, "\n");
    }
    fclose(pgmimg);
	
}