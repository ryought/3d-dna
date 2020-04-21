#### Description: Script to build raw contact count graph.
#### Written by: OD

# gawk -v P0=-2.148477541013386  -v P1=-0.9848534664777571 -v K0=3162 -v K=75000 -f test.awk test.data
function p(d)
{
  # p(d) =
  #     (d<K0) | PK0 - PK
  #     (d>K)  | PK  - PK
  # (otherwise)| P1 log(d) + P0 - PK
  return (d<K0) ? PK0PK : ((d>K) ? 0 : P0PK + (P1 * log(d)))
}
function abs(value)
{
	return (value<0?-value:value);
}
BEGIN {
  PK  = P0 + (P1 * log(K))
  PK0 = P0 + (P1 * log(K0))
  PK0PK = PK0 - PK
  P0PK = P0 - PK
}
# read in the cprops
{
	if (FILENAME==ARGV[1])
	{
		cname[$1]=$2
		clength[$2]=$3
		next
	}
}
# read in current assembly
{
	if (FILENAME==ARGV[2])
	{
		n=split($0,a)
		slength[FNR]=0
		for (k=1; k<=n; k++)
		{
			slist[abs(a[k])]=FNR
			shift[abs(a[k])]=slength[FNR]
			slength[FNR]+=clength[abs(a[k])]
			if (a[k]==abs(a[k]))
				orientation[a[k]]=0
			else
				orientation[abs(a[k])]=1
		}
		next
	}
}
{
	#check that both reads don't fall on the same contig and have been aligned with MAPQ>=MAPQ variable - optional
	if ((cname[$2] in slist) && (cname[$6] in slist) && (slist[cname[$2]]!=slist[cname[$6]]) && ($9 >= MAPQ) && ($12 >= MAPQ))
	{
		$2=cname[$2]
		$6=cname[$6]
		#recalculate positions
		if (orientation[$2] == 0)
			pos1 = $3 + shift[$2]
		else
			pos1 = -$3 + shift[$2] + clength[$2] + 1
		if (orientation[$6] == 0)
			pos2 = $7 + shift[$6]
		else
			pos2 = -$7 + shift[$6] + clength[$6] + 1

    # scaffold length
    slen1 = slength[slist[$2]]
    slen2 = slength[slist[$6]]

    # separation distance

		if (slist[$2] < slist[$6])
		{
      # <--->
      dist = pos1 + pos2
      pd = p(dist)
      if (pd > 0) {
        count[slist[$2]" "slist[$6]" "1] += 1
        prob[slist[$2]" "slist[$6]" "1] += pd
      }

      # <-<-
      dist = pos1 + slen2 - pos2
      pd = p(dist)
      if (pd > 0) {
        count[slist[$2]" "slist[$6]" "2] += 1
        prob[slist[$2]" "slist[$6]" "2] += pd
      }

      # ->->
      dist = slen1 - pos1 + pos2
      pd = p(dist)
      if (pd > 0) {
        count[slist[$2]" "slist[$6]" "3] += 1
        prob[slist[$2]" "slist[$6]" "3] += pd
      }

      # -><-
      dist = slen1 - pos1 + slen2 - pos2
      pd = p(dist)
      if (pd > 0) {
        count[slist[$2]" "slist[$6]" "4] += 1
        prob[slist[$2]" "slist[$6]" "4] += pd
      }
		}
		else
		{
      dist = pos2 + pos1
      pd = p(dist)
      if (pd > 0) {
        count[slist[$6]" "slist[$2]" "1] += 1
        prob[slist[$6]" "slist[$2]" "1] += pd
      }

      dist = pos2 + slen1 - pos1
      pd = p(dist)
      if (pd > 0) {
        count[slist[$6]" "slist[$2]" "2] += 1
        prob[slist[$6]" "slist[$2]" "2] += pd
      }

      dist = slen2 - pos2 + pos1
      pd = p(dist)
      if (pd > 0) {
        count[slist[$6]" "slist[$2]" "3] += 1
        prob[slist[$6]" "slist[$2]" "3] += pd
      }

      dist = slen2 - pos2 + slen1 - pos1
      pd = p(dist)
      if (pd > 0) {
        count[slist[$6]" "slist[$2]" "4] += 1
        prob[slist[$6]" "slist[$2]" "4] += pd
      }
		}
	}

}
END{
for (var in count)
	{
		print var, prob[var], count[var]
	}
}
