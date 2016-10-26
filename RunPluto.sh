scriptname="run_pluto.sh" 
rm -f ${scriptname}
echo "#!/bin/bash" >> ${scriptname}
echo "#PBS -N PLUTO" >> ${scriptname}
echo "#PBS -l walltime=48:00:00" >> ${scriptname}
echo >> ${scriptname}
echo "cd WASA-preselection/WASA-libs" >> ${scriptname}
echo "./run_plutto over He3 eta" >> ${scriptname}
echo "./run_plutto all He3 pi0" >> ${scriptname}
echo "./run_plutto all He3 pi0 pi0" >> ${scriptname}
echo "./run_plutto all He3 pi0 pi0 pi0" >> ${scriptname}
echo >> ${scriptname}
echo "rm -f $PWD/${scriptname}" >> ${scriptname}
chmod u+x ${scriptname}
qsub ${scriptname}
echo "${scriptname} generated and executed"
