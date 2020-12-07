#!/bin/tcsh

set WDir = `pwd`
set tests = `\ls -l | grep ^d | awk '{print $NF}'`
set logfile = "$WDir/logfile"
set failed = "0"
set passed = "0"
set vaisd = "../../VAISD/bin/compute_specden"

if ( -e $logfile ) then
    rm $logfile
endif

foreach test ( $tests )

  cd $test

  echo "Running test:  $test" >> $logfile

  if ( -e test_force.log || -e test_force.out ) then

    set freqfile = `\ls *freq*`
    set forcefile = `\ls *force*`
    $vaisd --gs $freqfile --es $forcefile -o test --figext png

  else

    set freqfile = `\ls *freq*`
    set esfile = `\ls *esopt*`
    $vaisd -m AS --gs $freqfile --es $esfile -o test --figext png

  endif

  diff test.results.txt ref.results.txt >& /dev/null
  
  if ( $status != "0" ) then
    echo "Test $test FAILED!" >> $logfile
    echo >> $logfile
    @ failed ++ 
  else
    echo "Test $test PASSED!" >> $logfile
    echo >> $logfile
    @ passed ++ 
  endif

  cd ..

end

# set gp = `which gnuplot`
# 
# if ( $? == "0" ) then
#   $gp --persist plot.gp
# endif
