void External(void)
{
  
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {   
      h1[point][0] += H_EXT_1;
      h2[point][0] += H_EXT_2;
      h3[point][0] += H_EXT_3;
    }
}
