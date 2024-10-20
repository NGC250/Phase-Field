void Laplacian(double* array, double* L_phi){
    
    double diff_area = 1/(step * step);
    double up, down, left, right, center;

    for (int i = 0; i < N; i++){
        
        up = array[(i >= H) ? (i - H) : (N - H + i)];
        down = array[(i + H < N) ? (i + H) : (i + H - N)];

        left = array[(i % W == 0) ? (i + W - 1) : (i - 1)];
        right = array[((i + 1) % W == 0) ? (i + 1 - W) : (i + 1)];

        center = array[i];

        L_phi[i] = ((right + left - 2 * center) + (up + down - 2 * center)) * diff_area;
    }
}