void initialise_tensors(void)
{
  for(p = 0; p < 3; p++)
    {
      for(q = 0; q < 3; q++)
        {
	  Delta[p][q] = (p == q) ? 1.0 : 0.0;

	  Mstr_coeff[p][q] = 1.5 * (l100 * Delta[p][q] + l111 * (1.0 - Delta[p][q]));
        }
    }

  for(p = 0; p < 3; p++)
    {
      for(q = 0; q < 3; q++)
        {
	  for(r = 0; r < 3; r++)
            {
	      for(s = 0; s < 3; s++)
                {                    
		  if((p == q) && (r == s))
                    {
                      if(q == r) Stiffness[p][q][r][s] = c11;
                      else Stiffness[p][q][r][s] = c12;
                    }
                    
		  else if((p != q) && (r != s))
                    {
		      if((p == r) && (q == s)) Stiffness[p][q][r][s] = c44;
                    }
                    
		  else Stiffness[p][q][r][s] = 0.0;
                }
            }
        }
    }
}
