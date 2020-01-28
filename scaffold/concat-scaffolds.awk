# concatenate all scaffolds
{
  if (all) {
    all=all" "$0
  } else {
    all=$0
  }
}
END {
  print all
}
