#!/usr/local/bin/vmd -dispdev text

# compute electrostatic potential map using the MEAD program package
# (Don Bashford, TSRI) and place positive and negative ions to
# potential minima and maxima, respectively, in random order. Maintain
# minimum ion-to-molecule and ion-to-ion distances.
# (C) Ilya Balabin <ilya.balabin@duke.edu>, Duke University, 2004-2006

### Tell Tcl that we're a package and any dependencies we may have
package require psfgen
package provide meadionize 1.2

# meadionize syntax
proc meadionize_usage {} {
    puts "Meadionize usage:"
    puts "meadionize -psf file.psf -pdb file.pdb -par par_all27_prot_na.inp -ipos na ineg cl -is 0.05"
    puts ""
    puts "Mandatory parameters:"
    puts "  -psf <value> - PSF file with solvated molecule"
    puts "  -pdb <value> - PDB file with solvated molecule"
    puts "  -par <value> - Charmm format parameter file"
    puts "  -ipos <value> - positive ion type (Na, K, Mg, Ca)"
    puts "  -ineg <value> - negative ion type (Cl)"
    puts "  Use ONLY ONE of the following two options:"
    puts "  -is <value> - required ionic strength (moles/liter), OR"
    puts "  -npos <value> -nneg <value> - \#\'s of positive and negative ions."
    puts ""
    puts "Optional parameters (do not change these unless you have a reason):"
    puts "  -o <value> - output PSF and PDB file prefix (default \"ionized\");"
    puts "  -seg <value> - ion segment name (default \"ION\");"
    puts "  -from <value> - min distance between ion and molecule, A (default 5);"
    puts "  -between <value> - min distance between different ions, A (default 5);"
    puts "  -epsin <value> - dielectric constant inside molecule (default 2);"
    puts "  -epsext <value> - dielectric constant of solvent (default 80);"
    puts "  -solrad <value> - solvent probe radius, A (default 1.4);"
    puts "  -sterln <value> - ion exclusion layer thickness, A (default 2);"
    puts "  -gr <value> - grid resolution, A (default 2.5);"
    puts "  -T <value> - temperature, K (default 300)."
    puts ""
    puts "On low memory machines, one may need to increase gr to 2-5."
    error ""
}

### put command line options in global array cmdline
proc cmdoptions { args } {
    global cmdline
    set args [eval "concat $args"]
    set n [llength $args]
    array unset cmdline
    for { set i 0 } { $i < $n } { incr i 2 } {
        set key [lindex $args $i]
        set val [lindex $args [expr $i + 1]]
        set cmdline($key) $val
    }
    return 0
}

### check command line for mandatory parameters
proc check_params { flags } {
    global cmdline
    foreach f $flags {
        if { ![info exists cmdline($f)] } { return 1 }
    }
    return 0
}

### find an executable (or exit if not found)
proc which { name } {
    set pathname ""
    global env
    set dirs [split $env(PATH) ":"]
    foreach dir $dirs {
        if { [file executable "${dir}/${name}"] } {
            set pathname "${dir}/${name}"
            break
        }
    }
    if { $pathname == "" } {
        puts "ERROR: \"$name\" not found, exiting"
        quit
    }
    return $pathname
}

# read atom radii from parameter file
proc get_radii { parfile } {
    set rel 0
    set fh [ open "$parfile" "r" ]
    foreach line [split [read $fh] "\n"] {
	# radii records begin after this line...
	if { [regexp {^\!\s*atom.*ignored.*epsilon.*Rmin\/2.*ignored.*eps.*$} $line] } {
	    set rel 1
	}
	# ... and end at NBFIX or HBON line or end of file
	if { [regexp {^\s*NBFIX.*$} $line] || [regexp {^\s*HBON.*$} $line] } {
	    set rel 0
	}
	# skip line if not relevant or comment
	if { !$rel || [regexp {^\s*\!.*$} $line] || [regexp {^\s*$} $line] } {continue}
	# parse line
	set fields {}
	foreach el [split $line] {
	    if {$el == {}} {continue}
	    lappend fields $el
	}
	set t [lindex $fields 0]
	set r [lindex $fields 3]
	# puts "  radius($t) = $r"
	set radii($t) $r
    }
    close $fh
    return [array get radii]   ;# array converted into list to return
}

# generate OGM file (grid dimension and resolution)
proc generate_ogm { topmol gr prefix } {
    # measure system dimensions
    set sel [atomselect $topmol all]
    foreach {min max} [measure minmax $sel] {break}
    $sel delete
    # find largest dimension
    set dim 0
    foreach r [vecsub $max $min] {
	if {$r > $dim} {set dim $r}
    }
    # number of grid nodes (with 5 extra nodes either side)
    set n [expr int([ expr $dim/$gr/2])*2 + 11]
    set fh [open "${prefix}.ogm" "w"]
    puts $fh "ON_ORIGIN ${n} ${gr}"
    # puts $fh "ON_CENT_OF_INTR ${n} ${gr}"
    close $fh
    return $dim
}

# generate PQR file (charges and radii, all atoms)
proc generate_pqr { topmol rlist pdbfile prefix } {

    # get atomic radii and charges; write temp PQR file
    array set radii $rlist   ;# convert list back to array
    set sel [atomselect $topmol "not water and not name SOD POT CAL MG CLA"]
    set reslist [$sel get resid]
    set charges [$sel get charge]
    set lradii {}
    foreach type [$sel get type] {
	lappend lradii $radii($type)
    }
    $sel writepdb "${prefix}.pqr"
    $sel delete

    # read ATOM/HETATM records from temp PQR file
    set fh [open "${prefix}.pqr" "r"]
    foreach line [split [read $fh] "\n"] {
	if { ![regexp {^ATOM\s+.*$} $line] } {continue}
	lappend pdblines $line
    }
    close $fh

    # consistency check
    if { [llength $pdblines] != [llength $reslist] } {
	error "INTERNAL ERROR: number of PDB/PQR lines different from number of atoms"
	quit
    }

    # set new fields and write real PQR file
    set fh [open "${prefix}.pqr" "w"]
    set n 0
    set resref 0
    set aind 1
    foreach line $pdblines res $reslist c $charges r $lradii {
	if { $resref != $res } { set resref $res ; incr n }
	set tstr " "; append tstr [format "%5d" $n]
	set line [string replace $line 21 26 $tstr]

	# rebuild index for structures over 99,999 atoms
	set line [string replace $line 4 10 [format "%7d" $aind]]
	incr aind

	set line [string range $line 0 53]
	append line [format " %6.3f" $c] [format " %6.3f" $r]
	puts $fh $line
    }

    close $fh
    return 0
}

# generate FPT file (coordinates of water oxygens that can be replaced)
proc generate_fpt { topmol from between prefix } {

    set ions "name SOD POT CAL MG CLA"
    set mol0 "(not water) and (not ${ions})"
    set selstr "noh and water"
    set selstr "${selstr} and not (within $from of (${mol0}))"
    set selstr "${selstr} and not (within $between of ${ions})"
    set sel [atomselect $topmol "${selstr}"]
    set watlist [$sel get index]
    set xyz [$sel get {x y z}]
    $sel delete
    set ffpt [open "${prefix}.fpt" "w"]

    foreach pos $xyz {
	foreach {x y z} $pos {break}
	set sx [format "%7.3f" $x]
	set sy [format "%7.3f" $y]
	set sz [format "%7.3f" $z]
	puts $ffpt "$sx  $sy  $sz"
    }

    close $ffpt
    return $watlist
}

# get ionic strength from numbers of positive and negative ions
proc nums2is { mol npos zpos nneg zneg } {

    # find number of water molecules
    set sel [atomselect $mol "water and noh"]
    set nWater [$sel num]
    puts "Found ${nWater} water molecules"
    $sel delete

    # compute ionic strength
    set nIS [expr ($npos*$zpos*$zpos + $nneg*$zneg*$zneg)]
    set is [expr $nIS / 0.01846 / $nWater]
    return $is

}

# get numbers of positive and negative ions from ionic strength
proc is2nums { mol is ipos zpos ineg zneg } {

    # compute initial net charge
    set sel [atomselect $mol all]
    set netCharge [eval "vecadd [$sel get charge]"]
    $sel delete
    puts "Initial net charge: ${netCharge}e"

    # find number of water molecules
    set sel [atomselect $mol "water and noh"]
    set nWater [$sel num]
    $sel delete
    puts "Found ${nWater} water molecules"

    # compute numbers of ions
    puts "Requested ion strength ${is} moles"
    set nIS [expr 0.01846 * $is * $nWater]
    set npos [expr (2.0 * $nIS - $netCharge*$zneg) / ($zpos * ($zpos + $zneg))]
    set npos [expr int($npos + 0.5)]
    set nneg [expr ($netCharge + $npos*$zpos) / $zneg]    
    set nneg [expr int($nneg + 0.5)]

    # check if ionic strength is sufficient for given net charge and water volume
    if { $npos <= 0 } {
	set npos 0
	set nneg [expr int($netCharge/$zneg + 0.5)]
	set is [nums2is $mol $npos $zpos $nneg $zneg]
	puts "WARNING: Ionic strength too low, cannot add $ipos ions!"
    }
    if { $nneg <= 0 } {
	set npos [expr int(-$netCharge/$zpos + 0.5)]
	set nneg [expr ($netCharge + $npos*$zpos) / $zneg]    
	set nneg [expr int($nneg + 0.5)]
	set is [nums2is $mol $npos $zpos $nneg $zneg]
	puts "WARNING: Ionic strength too low, cannot add $ineg ions!"
    }
    return [list $npos $nneg]

}

# main procedure
proc meadionize { args } {

    # internal
    set ph "MEADIONIZE)"

    # find MEAD executable
    set potential [which "potential"]
    puts "$ph Found MEAD utility $potential"

    # get command line options
    global cmdline
    cmdoptions $args

    # sanity check
    if { [check_params {"-psf" "-pdb" "-par" "-ipos" "-ineg"}] } {
        meadionize_usage; return 1
    }
    
    # set mandatory parameters
    set psffile $cmdline(-psf)
    set pdbfile $cmdline(-pdb)
    set parfile $cmdline(-par)
    if { [info exists cmdline(-is)] } {
	set is $cmdline(-is)
    } elseif { [info exists cmdline(-npos)] && [info exists cmdline(-nneg)] } {
	set npos $cmdline(-npos)
	set nneg $cmdline(-nneg)
    } else {
        meadionize_usage; return 1
    }
    
    # set ion types and charges (THESE MUST BE IN SYNC WITH THE ION TOPOLOGY FILE)
    set ipos $cmdline(-ipos)
    if { [regexp -nocase {^na$} $ipos] } {
	set ipos "SOD"
	set zpos 1
    } elseif { [regexp -nocase {^k$} $ipos] } {
	set ipos "POT"
	set zpos 1
    } elseif { [regexp -nocase {^ca$} $ipos] } {
	set ipos "CAL"
	set zpos 2
    } elseif { [regexp -nocase {^mg$} $ipos] } {
	set ipos "MG"
	set zpos 2
    } else {
	meadionize_usage; return 1
    }
    set ineg $cmdline(-ineg)
    if { [regexp -nocase {^cl$} $ineg] } {
	set ineg "CLA"
	set zneg 1
    } else {
	meadionize_usage; return 1
    }

    # set optinal parameters
    set outfile "ionized"; if {[info exists cmdline(-o)]} {set outfile $cmdline(-o)}
    set ionseg "ION"; if {[info exists cmdline(-seg)]} {set ionseg $cmdline(-seg)}
    set from 5; if {[info exists cmdline(-from)]} {set from $cmdline(-from)}
    set between 5; if {[info exists cmdline(-between)]} {set between $cmdline(-between)}
    set epsin 2.0; if {[info exists cmdline(-epsin)]} {set epsin $cmdline(-epsin)}
    set epsext 80.0; if {[info exists cmdline(-epsext)]} {set epsext $cmdline(-epsext)}
    set solrad 1.4; if {[info exists cmdline(-solrad)]} {set solrad $cmdline(-solrad)}
    set sterln 2.0; if {[info exists cmdline(-sterln)]} {set sterln $cmdline(-sterln)}
    set gr 2.5; if {[info exists cmdline(-gr)]} {set gr $cmdline(-gr)}
    set T 300; if {[info exists cmdline(-T)]} {set T $cmdline(-T)}

    # set package WD and files
    global env
    if ([info exists env(MEADIONIZEDIR)]) {
	set topfile $env(MEADIONIZEDIR)/ions.top
    } else {
	set topfile [file normalize [file dirname [info script]]]/ions.top
    }
    
    # print parameters
    puts "$ph Starting meadionize:"
    puts "$ph   PSF file ${psffile}"
    puts "$ph   PDB file ${pdbfile}"
    puts "$ph   Charmm ion topology file ${topfile}"
    puts "$ph   Charmm parameter file ${parfile}"
    puts "$ph   positive ions ${ipos}"
    puts "$ph   negative ions ${ineg}"
    puts "$ph   output file prefix ${outfile}"
    puts "$ph   added ions segment name ${ionseg}"
    puts "$ph   min distance between ion and molecule ${from} A"
    puts "$ph   min distance between different ions ${between} A"
    puts "$ph   dielectric constant inside molecule ${epsin}"
    puts "$ph   dielectric constant in solvent ${epsext}"
    puts "$ph   solvation radius ${solrad} A"
    puts "$ph   ionic exclusion layer thickness ${sterln} A"
    puts "$ph   grid resolution ${gr} A"
    puts "$ph   temperature ${T} K"

    # check if all necessary files exist
    foreach name [list $psffile $pdbfile $parfile ] {
	if { ![file exists $name] || ![file size $name] } {
	    error "ERROR: File ${name} not found or is empty"
	    meadionize_usage; return 1
	}
    }

    # read in topology (unless read already)
    global TOPOLOGY_READ
    if { [catch {if {$TOPOLOGY_READ != 1} {}}] } {
        set TOPOLOGY_READ 1
        topology $topfile
    }
    
#     # load and mass-center molecule
#     puts "$ph Loading and centering ${psffile}/${pdbfile}"
#     mol load psf "$psffile" pdb "$pdbfile"
#     set topmol [molinfo top]
#     set sel [atomselect $topmol all]
#     $sel moveby [vecinvert [measure center $sel]]
#     $sel writepdb "${outfile}.pdb"
#     $sel delete
#     mol delete all

    # load and _geometrically_ center molecule (not mass-center)
    puts "$ph Loading and centering ${psffile}/${pdbfile}"
    mol load psf "$psffile" pdb "$pdbfile"
    set topmol [molinfo top]
    set sel [atomselect $topmol all]
    foreach {min max} [measure minmax $sel] {break}
    $sel moveby [vecinvert [vecscale 0.5 [vecadd $min $max]]]
    $sel writepdb "${outfile}.pdb"
    $sel delete
    mol delete all

    # re-load centered molecule
    puts "$ph Re-loading ${psffile}/${pdbfile}"
    resetpsf
    mol delete all
    readpsf "$psffile"
    coordpdb "${outfile}.pdb"
    mol load psf "$psffile" pdb "${outfile}.pdb"
    set topmol [molinfo top]

    # check if ion segment already exists
    puts "$ph Checking if segment ${ionseg} already exists... "
    set sel [atomselect $topmol "segid $ionseg"]
    set n [llength [$sel get segid]]
    if { $n } {
	error "$ph yes; please use another ion segment name"
	quit
    } else {
	puts "$ph no, will add it"
    }
    $sel delete

    # compute npos and nneg from ionic strength or vice versa
    if { [info exists cmdline(-is)] } {
	foreach {npos nneg} [is2nums $topmol $is $ipos $zpos $ineg $zneg] {break}
    } else {
	set is [nums2is $topmol $npos $zpos $nneg $zneg]
    }
    set nIons [expr $npos + $nneg]    
    puts "$ph Will add $nIons ($npos $ipos and $nneg $ineg) ions to i.s. $is mol"

    # randomize order of adding ions (positive vs negative)
    puts "$ph Randomizing ion order"
    set templist {}
    for { set i 0 } { $i < $npos } { incr i } { lappend templist 1 }
    for { set i 0 } { $i < $nneg } { incr i } { lappend templist -1 }
    set ionlist {}
    while { [llength $templist] } {
	set thisNum [expr int([expr [llength $templist] * rand()])]
	lappend ionlist [lindex $templist $thisNum]
	set templist [lreplace $templist $thisNum $thisNum]
    }

    # get atom radii (converted back to array)
    puts "$ph Getting atom radii from ${parfile}"
    array set radii [get_radii $parfile]

    # generate OGM file
    puts "$ph Generating OGM file ${outfile}.ogm for MEAD"
    generate_ogm $topmol $gr $outfile

    # generate PQR file
    puts "$ph Generating PQR file ${outfile}.pqr for MEAD"
    generate_pqr $topmol [array get radii] $pdbfile $outfile
    
    # generate FPT file
    puts "$ph Generating FPT file ${outfile}.fpt for MEAD"
    set watlist [generate_fpt $topmol $from $between $outfile]
    
    # run MEAD
    puts "$ph Computing potential map; this may take a while..."
    set comlst [list "$potential" "-epsin" "$epsin" "-ionicstr" "$is"]
    if { [info exists cmdline(-epsext)] } {lappend comlst "-epsext" "$epsext"}
    if { [info exists cmdline(-solrad)] } {lappend comlst "-solrad" "$solrad"}
    if { [info exists cmdline(-sterln)] } {lappend comlst "-sterln" "$sterln"}
    if { [info exists cmdline(-T)] } {lappend comlst "-T" "${T}"}
    lappend comlst "$outfile" ">" "${outfile}.out"
    puts "$ph Running \"$comlst\""
    if { [catch {eval "exec $comlst"} pout ] } {
	puts "$ph WARNING: MEAD potential generated warnings/errors:\n$pout"
    }
    
    # associate potentials with water oxygen indices
    puts "$ph Potential map computed, reading potentials"
    set potlist {}
    set fh [open "${outfile}.out" "r"]
    set ignore 1
    foreach line [split [read $fh] "\n"] {
	# check whether reached actual potentials
	if { [regexp {^Only one level:.*$} $line] } {
	    set ignore 0 ; continue
	}
	if {$ignore} {continue}
	# check whether reached EOF
	if { [scan $line "%g" p] != 1 } { break }
	lappend potlist $p
    }
    close $fh

    # consistency check
    puts "$ph Associating potentials with water molecules"
    if { [llength $potlist] != [llength $watlist] } {
	error "$ph ERROR: number of potentials different from number of waters!"
	quit
    }
    foreach i $watlist p $potlist {set potmap($i) $p}


    # find water molecules to replace with ions
    puts "$ph Finding water molecules to replace with ions"
    set wat2pos {}
    set wat2neg {}

    for { set i 0 } { $i < $nIons } { incr i } {

	# find water to be replaced
	set mol0 "(not water) and (not name SOD POT CAL MG CLA)"
	set selstr "(noh and water)"
	set selstr "${selstr} and (not within ${from} of (${mol0}))"
	set selstr "${selstr} and not (within $between of name SOD POT CAL MG CLA)"
	if { [llength $wat2pos] } {
	    set selstr "${selstr} and (not within $between of (index ${wat2pos}))"
	}
	if { [llength $wat2neg] } {
	    set selstr "${selstr} and (not within $between of (index ${wat2neg}))"
	}
	
	set sel [atomselect $topmol "${selstr}"]
	set watlist [$sel get index]
	$sel delete

	# find potential min or max
	if { [lindex $ionlist $i] > 0 } {

	    # find potential min
	    set num [lindex $watlist 0]
	    set pot $potmap($num)
	    foreach j $watlist {
		if { $potmap($j) < $pot } {
		    set num $j
		    set pot $potmap($j)
		}
	    }
	    set sel [atomselect $topmol "index $num"]
	    set seg [$sel get segid]
	    set res [$sel get resid]
	    set name [$sel get name]
	    foreach {x y z} [lindex [$sel get {x y z}] 0] {break}
	    set x [format "%7.3f" $x]
	    set y [format "%7.3f" $y]
	    set z [format "%7.3f" $z]
	    set pot [format "%6.3f" $pot]
	    puts "  $num ${seg}:${res}:${name} at ${x} ${y} ${z} (${pot} eV) ==> $ipos"
	    $sel delete
	    lappend wat2pos $num

	} else {

	    # find potential max
	    set num [lindex $watlist 0]
	    set pot $potmap($num)
	    foreach j $watlist {
		if { $potmap($j) > $pot } {
		    set num $j
		    set pot $potmap($j)
		}
	    }
	    set sel [atomselect $topmol "index $num"]
	    set seg [$sel get segid]
	    set res [$sel get resid]
	    set name [$sel get name]
	    foreach {x y z} [lindex [$sel get {x y z}] 0] {break}
	    set x [format "%7.3f" $x]
	    set y [format "%7.3f" $y]
	    set z [format "%7.3f" $z]
	    set pot [format "%6.3f" $pot]
	    puts "  $num ${seg}:${res}:${name} at ${x} ${y} ${z} (${pot} eV) ==> $ineg"
	    $sel delete
	    lappend wat2neg $num
	}
    }

    # consistency check
    set posnum [llength $wat2pos]
    set negnum [llength $wat2neg]
    puts "$ph Will replace ${posnum} waters with ${ipos} and ${negnum} waters with ${ineg}"

    # getting water positions and deleting them
    puts "$ph Deleting water molecules to be replaced"
    if { $posnum } {
	set sel [atomselect $topmol "index $wat2pos"]
	set posposlist [$sel get { x y z}]
	set seglist [$sel get segid]
	set reslist [$sel get resid]
	foreach seg $seglist res $reslist {delatom $seg $res}
	$sel delete
    }
    if { $negnum } {
	set sel [atomselect $topmol "index $wat2neg"]
	set negposlist [$sel get { x y z}]
	set seglist [$sel get segid]
	set reslist [$sel get resid]
	foreach seg $seglist res $reslist {delatom $seg $res}
	$sel delete
    }

    # replacing selected water molecules with ions
    puts "$ph Replacing deleted water molecules with ions"
    set posreslist {}
    set negreslist {}
    for { set i 1 } { $i <= $npos } { incr i } {lappend posreslist $i}
    for { set i [expr $npos+1] } { $i <= $nIons } { incr i } {lappend negreslist $i}
    segment "$ionseg" {
	first NONE
	last NONE
	foreach i $posreslist {residue $i "${ipos}"}
	foreach i $negreslist {residue $i "${ineg}"}
    }

    # assign coordinates
    puts "$ph Assigning ion coordinates"
    if { $posnum } {
	foreach res $posreslist pos $posposlist {
	    coord $ionseg $res "${ipos}" $pos
	}
    }
    if { $negnum } {
	foreach res $negreslist pos $negposlist {
	    coord $ionseg $res "${ineg}" $pos
	}
    }
    
    # save new structure
    puts "$ph saving ${outfile}.psf/${outfile}.pdb"
    writepsf "${outfile}.psf"
    writepdb "${outfile}.pdb"
    
    puts "$ph All done."
    return 0

}

