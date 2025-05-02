typedef struct
{
    uint64_t max_index;
    uint64_t min_index;
    double max_value;
    double min_value;
}
Extract;

Extract In_Array(double* array, uint64_t size)
{
    Extract extremes;

    extremes.min_value = array[0];
    extremes.max_value = array[0];
    extremes.min_index = 0;
    extremes.max_index = 0;

    for(i = 1; i < size; i++)
    {
        if(array[i] < extremes.min_value)
        {
            extremes.min_value = array[i];
            extremes.min_index = i;
        }
        if(array[i] > extremes.max_value)
        {
            extremes.max_value = array[i];
            extremes.max_index = i;
        }
    }
    return extremes;
}

uint64_t* Voronoi(uint64_t num_cells)
{
    uint64_t elements , rand_x , rand_y;
    uint64_t *tessellation , cell_centers[2 * num_cells];
    double x , y , x_c , y_c , distances[num_cells];

    if((tessellation = (uint64_t*)malloc(N * sizeof(uint64_t))) == NULL)
    {
        printf("Tessellation could not be allocated!"); exit(1);
    }

    elements = 0;
    while(elements < num_cells)
    {
        rand_x = rand() % H;
        rand_y = rand() % W;
        cell_centers[elements] = rand_x;
        cell_centers[elements + num_cells] = rand_y;

        elements++;
    }

    Extract extremes;

    for(uint64_t i = 0; i < H; i++)
    {
        for(uint64_t j = 0; j < W; j++)
        {
            x = (double)(i);
            y = (double)(j);
            
            for(grain = 0; grain < num_cells; grain++)
            {
                x_c = (double)(cell_centers[grain]);
                y_c = (double)(cell_centers[grain + num_cells]);
                
                distances[grain] = (x - x_c) * (x - x_c) + (y - y_c) * (y - y_c);
            }
            
            extremes = In_Array(distances, num_cells);
            
            tessellation[i * W + j] = extremes.min_index;
        }
    }

    return tessellation;
}
