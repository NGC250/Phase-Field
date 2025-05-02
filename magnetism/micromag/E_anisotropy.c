void Anisotropy()
{
  double m1_sq , m2_sq , m3_sq;
  
  for(point = 0; point < N; point++)
    {
      m1_sq = m1[point][0] * m1[point][0];
      m2_sq = m2[point][0] * m2[point][0];
      m3_sq = m3[point][0] * m3[point][0];

      h1[point][0] -= m1[point][0] * (K1_S * (1.0 - m1_sq) + K2_S * m2_sq * m3_sq);
      h2[point][0] -= m2[point][0] * (K1_S * (1.0 - m2_sq) + K2_S * m3_sq * m1_sq);
      h3[point][0] -= m3[point][0] * (K1_S * (1.0 - m3_sq) + K2_S * m1_sq * m2_sq);
    }
}
