FUNCTION doslopefunc,num,x,y
  avex = 0.0
  avey = 0.0
  sum1 = 0.0
  sum2 = 0.0
  for i=0,num-1 do begin
    avex += x[i]
    avey += y[i]
  endfor
  avex /= num
  avey /= num
  for i=0,num-1 do begin
    sum1 += (x[i]-avex)*(y[i]-avey)
    sum2 += (x[i]-avex)*(x[i]-avex)
  endfor
  return, sum1/sum2
END
