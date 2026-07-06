# awk script to postprocess the data measurements

# call like:
# cat results-all-reformatted.output | gawk -M -f datapostproc.awk | sort -n -k 2 -k 1 > results-all-avg-min-max.data
#
# we have a data file with several values for every measurement
# the first column contains a key per measurement such as:
# method1-xvalue1 10 1.5
# method1-xvalue1 10 1.6
# ...
# method1-xvalue2 20 1.2
# method1-xvalue2 20 1.3
#
# the ordering of the lines is irrelevant and can be mixed

function myminf(arr, key, ind) {
  # minimum
  # per new key, take the first value. otherwise, take the minimum
  if (!($key in arr) || $ind < arr[$key])
    arr[$key] = $ind
}

function mymaxf(arr, key, ind) {
  # maximum
  # per new key, take the first value. otherwise, take the maximum
  if (!($key in arr) || $ind > arr[$key])
    arr[$key] = $ind
}

function mysumf(arr, key, ind) {
  # sum up
  arr[$key] += $ind
}

BEGIN {
  PREC=100
}

# start per input line processing
{
  # define what is where in the source data file
  key                    = 1
  size                   = 2
  iterationsuntiltimeout = 3
  nsperop                = 4
  nsperopunit            = 5
  speed                  = 6
  speedunit              = 7
  arch                   = 8
  archvalue              = 9
  os                     = 10
  osvalue                = 11
  threads                = 12
  threadsvalue           = 13

  if (!($key in cnt)) {
    orderindex++;
    order[orderindex] = $key;
  }

  # count the observed repetitions of each key
  cnt[$key]++;

  # keep last values of some fields
  arr_size[$key] = $size
  arr_speedunit[$key] = $speedunit
  arr_arch[$key] = $arch
  arr_archvalue[$key] = $archvalue
  arr_os[$key] = $os
  arr_osvalue[$key] = $osvalue
  arr_threads[$key] = $threads
  arr_threadsvalue[$key] = $threadsvalue

  # minimum
  # per new key, take the first value. otherwise, take the minimum
  # if (!($1 in mymin) || $3 < mymin[$1])
  #   mymin[$1] = $3

  myminf(mymin, key, iterationsuntiltimeout)
  myminf(myspeedmin, key, speed)

  # maximum
  # if (!($1 in mymax) || $3 > mymax[$1])
  #   mymax[$1] = $3
  mymaxf(mymax, key, iterationsuntiltimeout)
  mymaxf(myspeedmax, key, speed)

  # sum (for mean)
  mysumf(myspeedsum, key, speed)

  # median, percentiles
  # would need to remember all data points
  #  n = length(a)
  #  asort(a)
  #
  #  # Calculate median
  #  if (n % 2 == 0) {
  #      median = (a[n/2] + a[n/2 + 1]) / 2
  #  } else {
  #      median = a[(n+1)/2]
  #  }
  #
  #  # Calculate percentiles
  #  p25 = a[int(n*0.25)+1]
  #  p75 = a[int(n*0.75)+1]
  # #  # Print results
  #  print "Median:", median
  #  print "25th Percentile:", p25
  #  print "75th Percentile:", p75
}

END {
  # mean
  for (i in myspeedsum)
    myspeedmean[i] = myspeedsum[i]/cnt[i]

  # emit a header block:
  print "# column 01: key"
  print "# column 02: array length in byteskey"
  print "# column 03: 'speed:'"
  print "# column 04: arithmetical average speed"
  print "# column 05: minimum speed"
  print "# column 06: maximum speed"
  print "# column 07: 'MB/s' the unit of speed"
  print "# column 08: repetitions, how many data points where used for statistics"
  print "# column 09: 'arch:'"
  print "# column 10: processor architecture"
  print "# column 11: 'os:'"
  print "# column 12: operating system"
  print "# column 13: 'threads:'"
  print "# column 14: number of threads for parallel alg. all other use 1 thread"
  for (i=1; i <= length(order); i++) {
    key = order[i]
    print key, arr_size[key],

      ################################## speed
      "speed: " myspeedmean[key],
      # min, max
      myspeedmin[key], myspeedmax[key],
      arr_speedunit[key],

      cnt[key],
      # remaining fields from input
      arr_arch[key], arr_archvalue[key],
      arr_os[key], arr_osvalue[key],
      arr_threads[key], arr_threadsvalue[key]
  }
}


