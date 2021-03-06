#### Description: Script to build raw contact count graph.
#### Written by: OD

# gawk -v P0=-2.148477541013386  -v P1=-0.9848534664777571 -v K0=3162 -v K=75000 -f test.awk test.data
function p_slow(d)
{
  if (d < K0) {
    return P0 + (P1 * log(K0));
  } else if (d > K) {
    return P0 + (P1 * log(K));
  } else {
    return P0 + (P1 * log(d));
  }
}
function p(d)
{
  return (d<K0) ? PK0 : ((d>K) ? PK : P0 + (P1 * log(d)))
}
function abs(value)
{
	return (value<0?-value:value);
}
BEGIN {
  PK  = P0 + (P1 * log(K))
  PK0 = P0 + (P1 * log(K0))
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
      count[slist[$2]" "slist[$6]" "1] += 1
      dist = pos1 + pos2
      prob[slist[$2]" "slist[$6]" "1] += p(dist)

      # <-<-
      count[slist[$2]" "slist[$6]" "2] += 1
      dist = pos1 + slen2 - pos2
      prob[slist[$2]" "slist[$6]" "2] += p(dist)

      # ->->
      count[slist[$2]" "slist[$6]" "3] += 1
      dist = slen1 - pos1 + pos2
      prob[slist[$2]" "slist[$6]" "3] += p(dist)

      # -><-
      count[slist[$2]" "slist[$6]" "4] += 1
      dist = slen1 - pos1 + slen2 - pos2
      prob[slist[$2]" "slist[$6]" "4] += p(dist)
		}
		else
		{
      count[slist[$6]" "slist[$2]" "1] += 1
      dist = pos2 + pos1
      prob[slist[$6]" "slist[$2]" "1] += p(dist)

      count[slist[$6]" "slist[$2]" "2] += 1
      dist = pos2 + slen1 - pos1
      prob[slist[$6]" "slist[$2]" "2] += p(dist)

      count[slist[$6]" "slist[$2]" "3] += 1
      dist = slen2 - pos2 + pos1
      prob[slist[$6]" "slist[$2]" "3] += p(dist)

      count[slist[$6]" "slist[$2]" "4] += 1
      dist = slen2 - pos2 + slen1 - pos1
      prob[slist[$6]" "slist[$2]" "4] += p(dist)
		}
	}

}
END{
for (var in count)
	{
		print var, prob[var], count[var]
	}
}
