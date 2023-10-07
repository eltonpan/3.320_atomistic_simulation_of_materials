## Fixes for bugs for running `pw.x < pwscf.in > pwscf.out`

Problem 1: **pw.x command not found**
Fix: make sure to add executable directory to $PATH
`export PATH="/home/gridsan/epan1/3320_atomistic_shared/qe-7.2/bin:$PATH"`

Problem 2: Ge potential not available in shared folder. `wget https://www.quantum-espresso.org/upf_files/Ge.pz-dn-rrkjus_psl.0.2.2.UPF` doesn't work as the file is no longer available on the website:
Fix: Download Ge potential from another source: 
`cd ~/labutil/lab2_samples`
`wget https://pseudopotentials.quantum-espresso.org/upf_files/Ge.pz-dn-rrkjus_psl.0.2.2.UPF`
then change `pseudo_dir = '/home/gridsan/<USER>/3320_atomistic_shared/q-e/pseudo/'`in `pwscf.in` to pseudo_dir = '/home/gridsan/<USER>/labutil/lab2_samples'


Now running `pw.x < pwscf.in > pwscf.out` should work