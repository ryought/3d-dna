# concatenate all scaffolds (composed of >1 contigs)
{
  if (NF > 1) {
    # append to A
    if (A) {
      A=A" "$0
    } else {
      A=$0
    }
  } else {
    # append to B
    if (B) {
      B=B"\n"$0
    } else {
      B=$0
    }
  }
}
END {
  print A
  print B
}
