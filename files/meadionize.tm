<TeXmacs|1.0.6>

<style|article>

<\body>
  <with|font-base-size|14|<paragraph|Meadionize Plugin for VMD, version
  1.1.>>

  <paragraph|Why Meadionize?>Meadionize is an improved version of the
  <hlink|Autoionize|http://www.ks.uiuc.edu/Research/vmd/plugins/autoionize/>
  plugin for VMD. Autoionize randomly places sodium and chlorine counterions
  around a solvated molecule. Autoionize has been designed for molecules with
  small net charges, where the electrostatic interactions between the
  molecule and the counterions are relatively weak, and the molecular
  structure and function are not sensitive to the counterion distribution.
  That is not the case for highly charged systems, most notably nucleic
  acids, which are surrounded by a clowd of counterions; correct simulations
  of these systems require placing ions according to the electrostatic
  potential of the molecule. Meadionize addresses this problem by placing the
  ions into the minima of the electrostatic potential map generated by the
  potential utility of the <hlink|MEAD|ftp://ftp.scripps.edu/pub/bashford/>
  program package by <hlink|Don Bashford|http://www.scripps.edu/mb/bashford/>
  by solving the Poisson-Boltzmann equation. Meadionize accepts all MEAD
  configuration parameters, generates all necessary input files for the
  potential program, executes it, and uses its output to place the
  counterions.

  <paragraph|Installation.>Create a directory for Meadionize: e.g., under
  your home directory<next-line>

  <\with|font-family|tt>
    mkdir $HOME/meadionize
  </with>

  \;

  <no-indent>Download and save files <hlink|meadionize.tcl|http://www.chem.duke.edu/~ilya/software/Meadionize/meadionize.tcl>
  (the main script), <hlink|ions.top|http://www.chem.duke.edu/~ilya/software/Meadionize/ions.top>
  (the ion topology), and <hlink|pkgIndex.tcl|http://www.chem.duke.edu/~ilya/software/Meadionize/pkgIndex.tcl>
  (the package index that tells VMD where to look for the package) into that
  directory. Add the following two lines to your
  <with|font-family|tt|$HOME/.vmdrc> file (create one if you do not have it):

  \;

  <\with|font-family|tt>
    global env

    lappend auto_path $env(HOME)/meadionize
  </with>

  \;

  <no-indent>Then, download and install the
  <hlink|MEAD|ftp://ftp.scripps.edu/pub/bashford/> program package following
  the included installation instructions.

  <paragraph|Usage.>Like all other VMD plugins, you first need to load
  Meadionize into VMD. Type the following command in the VMD console (better
  yet, <hlink|Tk console|http://www.ks.uiuc.edu/Research/vmd/plugins/vmdtkcon/>):<line-break>

  \;

  <with|font-family|tt|package require meadionize>

  \;

  <no-indent>This command actually loads Meadionize into VMD and displays the
  installed plugin version. Run '<with|font-family|tt|meadionize>' with no
  parameters to get brief help on Meadionize syntax. Meadionize requires the
  following mandatory parameters: the PSF (structure) file name, the PDB
  (coordinate) file name, the Charmm parameter file name, positive and
  negative ion types (as of this point, <with|mode|math|Na<rsup|+>>,
  <with|mode|math|K<rsup|+>>, <with|mode|math|Ca<rsup|2+>>,
  <with|mode|math|Mg<rsup|2+>>, and <with|mode|math|Cl<rsup|->> ions are
  supported), and either ionic strength, or explicit numbers of the positive
  and negative ions to add. Other parameters are optional; their default
  values rarely need to be changed, except the <with|font-family|tt|fg>
  option (fine grid resolution), which may need to be increased to
  <with|mode|math|0.5 - 1.0><specific|latex|\\AA>A when running Meadionize on
  machines with low RAM and/or adding ions to an exceptionally large system.
  This is an example of Meadionize command line:<line-break>\ 

  \;

  <with|font-family|tt|meadionize -psf solvated.psf -pdb solvated.pdb -par
  par_all27_prot_lipid.prm -ipos na -ineg cl -is 0.1>

  \;

  <no-indent>Since calculating the electrostatic potential map for large
  molecules can take hours, it is recommended to run Meadionize using VMD in
  the text-only mode. To do that, one needs to copy the above two command
  lines into a file (e.g., <with|font-family|tt|do_ionize.tcl>), add a
  command <with|font-family|tt|quit> to the end of the file to tell VMD to
  stop after executing Meadionize, and run the following command in the UNIX
  shell (csh or tcsh):<line-break>\ 

  \;

  <with|font-family|tt|vmd -dispdev text \<less\> do_ionize.tcl \<gtr\>&
  do_ionize.log &>

  \;

  <no-indent>For bash or ksh, the shell command syntax is slightly
  different:<line-break>\ 

  \;

  <with|font-family|tt|vmd -dispdev text \<less\> do_ionize.tcl \<gtr\>
  do_ionize.log 2\<gtr\>&1 &>

  \;

  <no-indent>To monitor Meadionize in real time, use the following command
  (any shell):<line-break>\ 

  \;

  <with|font-family|tt|tail -f do_ionize.log>

  \;

  <paragraph|Output.>Meadionize prints out diagnostic messages about
  performed steps as well as repeats messages from the potential utility.
  Often, these messages include the following warnings:

  \;

  <\with|font-family|tt>
    WARNING: SAVanal_calc: vertex found with count = 2

    WARNING: SAVanal_calc: vertex found with count = 2

    WARNING: SAVanal_calc: vertex found with count = 1
  </with>

  \;

  <no-indent>These warnings are harmless, and arise from numerical
  degeneracies in the calculation of the molecular surface (see You and
  Bashford <with|font-shape|italic|J. Comp. Chem.>
  <with|font-series|bold|16>, 743 (1995)). Don Bashford explains these
  warnings in more details in his post on the <hlink|Computational Chemistry
  List|http://ccl.osc.edu/cgi-bin/ccl/message.cgi?2003+02+27+002>.

  <paragraph|Algorithm.>Meadionize performs the following basic steps:

  <\enumerate-numeric>
    <item>If an ionic strength is requested, Meadionize finds the number of
    water molecules and calculates the numbers of the positive and negative
    ions (alternatively, these numbers are given as parameters).

    <item>Prepares necessary input files and calls the
    <with|font-family|tt|potential> utility of the MEAD package to solve the
    Poisson-Boltzmann equation and calculate the electrostatic potential map.

    <item>In random order, replaces the water molecules at the electrostatic
    potential minima (for positive ions) or maxima (for negative ions) with
    the corresponding ions. Each time, the new ion is placed so that a
    minimum distance between any two ions, as well as a minimum distance
    between any ion and themolecule, are maintained.
  </enumerate-numeric>

  At the first step, Meadionize computes the numbers of positive and negative
  ions from two conditions: zero net charge of the system and the ionic
  strength:

  <\itemize>
    <item><with|mode|math|n<rsub|-> z<rsub|-> - n<rsub|+> z<rsub|+> =
    Q<rsub|mol>> , where <with|mode|math|n<rsub|->>
    (<with|mode|math|n<rsub|+>>) and <with|mode|math|z<rsub|->>
    (<with|mode|math|z<rsub|+>>) are numbers and charges of the negative
    (positive) ions, respectively, and <with|mode|math|Q<rsub|mol>> is the
    net charge of the molecule before adding ions.

    <item><with|mode|math|n<rsub|-> z<rsup|<rsup|>2><rsub|-> \ + n<rsub|+>
    z<rsup|2><rsub|+> = 2 N<rsub|is>> , where <with|mode|math|N<rsub|is>> is
    a quantity proportional to the total number of ions, a function of the
    requested ionic strength.\ 
  </itemize>

  In the second condition, <with|mode|math|N<rsub|is> = C \ N<rsub|A>
  \ V<rsub|H<rsub|2>O> > , where <with|mode|math|C> is the ionic strength
  (mol), <with|mode|math|N<rsub|A> = 6.022 \<cdot\> 10<rsup|23>> is the
  Avogadro number, and <with|mode|math|V<rsub|H<rsub|2>O>> is the water
  volume (L). The latter can be represented as
  <with|mode|math|><with|mode|math|V<rsub|H<rsub|2>O> = N<rsub|H<rsub|2>O> /
  \<rho\> = N<rsub|H<rsub|2>O > m<rsub|H<rsub|2>O> / \<rho\><rsub|m>> , where
  <with|mode|math|N<rsub|H<rsub|2>O>> is the number of water molecules in the
  system, <with|mode|math|\<rho\>> is the volume water density,
  <with|mode|math|\<rho\><rsub|m> = 0.982g/cm<rsup|3>> is the mass water
  density (TIP3 water model, Jorgensen <with|font-shape|italic|et al>,
  <with|font-shape|italic|J. Chem. Phys.> <with|font-series|bold|79>, 926
  (1983)), and <with|mode|math|m<rsub|H<rsub|2>O> = 18 \<cdot\> 1.67262
  \<cdot\> 10<rsup|-24>g> is the water molecule mass. Combining the above,
  one obtains the following formulas:\ 

  <\itemize>
    <item><with|mode|math|N<rsub|is> = 0.01846 \ \ C \ N<rsub|H<rsub|2>O>>

    <item><with|mode|math|n<rsub|+> = <mid|[>2 N<rsub|is> - Q<rsub|mol>
    z<rsub|-><right|]> <mid|/> <mid|[>z<rsub|+> <left|(>z<rsub|+> +
    z<rsub|-><right|)><mid|]>>

    <item><with|mode|math|n<rsub|-> = <left|[>Q<rsub|mol> + n<rsub|+>
    z<rsub|+><mid|]> <mid|/> z<rsub|->>
  </itemize>

  <paragraph|License.>Meadionize is released under the <hlink|GNU public
  license|http://www.gnu.org/copyleft/gpl.html>. I thank Don Bashford for
  developing MEAD, which made Meadionize possible. All comments, suggestions,
  and bug reports are very welcome.

  <no-indent>______________________________________________________________________

  <no-indent>(C) <hlink|Ilya Balabin|http://www.chem.duke.edu/~ilya>, Duke
  University, 2004-2006

  \;
</body>

<\initial>
  <\collection>
    <associate|font|helvetica>
    <associate|font-base-size|11>
    <associate|math-font|roman>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|1>>
    <associate|auto-4|<tuple|4|1>>
    <associate|auto-5|<tuple|5|2>>
    <associate|auto-6|<tuple|6|2>>
    <associate|auto-7|<tuple|7|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <with|par-left|<quote|6fn>|Meadionize Plugin for VMD, version 1.1.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Why Meadionize?
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Installation.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Usage. <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Output. <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Algorithm.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|License.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.15fn>>
    </associate>
  </collection>
</auxiliary>