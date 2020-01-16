# find minimum
{
  if (FNR==NR)
  {
    if (!min) {
      min = $4
    }
    if (min > $4) {
      min = $4
    }
  }
}
# output by setting offset
{
  if (FNR!=NR) {
    print $1, $2, $3, $4-min, $5
  }
}
