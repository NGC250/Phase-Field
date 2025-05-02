void savedata(uint32_t at)
{
  FILE *image_m;
  char filename[256];
	
  sprintf(filename , "../Data/status/domains%u.dat",at);
	
  image_m = fopen(filename , "w");

  fprintf(image_m , "mx , my , mz\n");
  for(point = 0; point < N; point++)
    {
      fprintf(image_m , "%le,%le,%le\n" , m1[point][0] , m2[point][0] , m3[point][0]);
    }
	
  fclose(image_m);
}
