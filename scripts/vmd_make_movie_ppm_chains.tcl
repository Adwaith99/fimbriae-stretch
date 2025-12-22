# Usage:
#   vmd -dispdev opengl -e vmd_make_movie_ppm_chains.tcl \
#     -args start.gro traj.xtc outprefix stride fps width height [chain_ranges.txt]

proc guessChainsFromGro {grofile} {
    # Heuristic: new "chain" whenever residue number decreases (resid reset)
    # Returns list of {start end} atom-index ranges (0-based, inclusive)
    set f [open $grofile r]
    gets $f title
    gets $f natomsLine
    set natoms [expr {[string trim $natomsLine]}]

    set prevResid ""
    set chainStarts [list 0]
    set atomIndex 0

    while {[gets $f line] >= 0} {
        if {$atomIndex >= $natoms} { break }

        # GRO resid is columns 1-5
        set residStr [string trim [string range $line 0 4]]
        if {$residStr ne ""} {
            set resid [expr {$residStr}]
            if {$prevResid ne "" && $resid < $prevResid} {
                lappend chainStarts $atomIndex
            }
            set prevResid $resid
        }
        incr atomIndex
    }
    close $f

    lappend chainStarts $natoms

    set ranges {}
    for {set i 0} {$i < [expr {[llength $chainStarts]-1}]} {incr i} {
        set s [lindex $chainStarts $i]
        set e [expr {[lindex $chainStarts [expr {$i+1}]] - 1}]
        lappend ranges [list $s $e]
    }
    return $ranges
}

proc readChainRangesFile {path} {
    # File format: each non-empty, non-# line has: start end
    set ranges {}
    set f [open $path r]
    while {[gets $f line] >= 0} {
        set line [string trim $line]
        if {$line eq ""} { continue }
        if {[string match "#*" $line]} { continue }
        set parts [split $line]
        if {[llength $parts] < 2} { continue }
        set s [expr {[lindex $parts 0]}]
        set e [expr {[lindex $parts 1]}]
        lappend ranges [list $s $e]
    }
    close $f
    if {[llength $ranges] == 0} {
        error "No ranges read from $path"
    }
    return $ranges
}

# -------- Args --------
set gro     [lindex $argv 0]
set xtc     [lindex $argv 1]
set outpre  [lindex $argv 2]
set stride  [expr {[lindex $argv 3]}]
set fps     [expr {[lindex $argv 4]}]
set width   [expr {[lindex $argv 5]}]
set height  [expr {[lindex $argv 6]}]
set rangesFile [lindex $argv 7]

# -------- Load --------
mol new $gro type gro waitfor all
set molid [molinfo top]
mol addfile $xtc type xtc waitfor all

# -------- Display settings --------
color Display Background white
axes location Off
display projection Orthographic
display resetview
display resize $width $height
display update

# Remove default rep
catch { mol delrep 0 $molid }

# -------- Chain ranges --------
if {$rangesFile ne "" && [file exists $rangesFile]} {
    set chainRanges [readChainRangesFile $rangesFile]
    puts "Info) Using manual chain ranges from: $rangesFile"
} else {
    set chainRanges [guessChainsFromGro $gro]
    puts "Info) Using auto chain ranges from GRO resid resets (n=[llength $chainRanges])"
}

# -------- Representations (one per chain chunk) --------
set colorCycle {0 1 3 4 5 6 7 8 9 10 11 12 13 14 15}
set rep 0
foreach r $chainRanges {
    lassign $r s e
    mol addrep $molid
    mol modselect $rep $molid "protein and index $s to $e"
    mol modstyle  $rep $molid NewCartoon 0.3 10.0 4.1 0
    set c [lindex $colorCycle [expr {$rep % [llength $colorCycle]}]]
    mol modcolor  $rep $molid ColorID $c
    incr rep
}

# -------- Render frames (PPM) --------
set nframes [molinfo $molid get numframes]
set framecount 0
for {set i 0} {$i < $nframes} {incr i $stride} {
    animate goto $i
    set img [format "%s_%05d.ppm" $outpre $framecount]
    render snapshot $img
    if {![file exists $img] || [file size $img] < 1000} {
        puts stderr "ERROR: frame not written or too small: $img"
        exit 1
    }
    incr framecount
}

# -------- Encode movie --------
set mp4 "${outpre}.mp4"
if {[catch {
    exec ffmpeg -y -framerate $fps -i "${outpre}_%05d.ppm" \
      -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
      -c:v libx264 -pix_fmt yuv420p $mp4
} err]} {
    puts stderr "ERROR: ffmpeg failed: $err"
    exit 1
}

if {![file exists $mp4] || [file size $mp4] < 10000} {
    puts stderr "ERROR: mp4 missing/too small after ffmpeg: $mp4"
    exit 1
}

puts "DONE: wrote $mp4"
quit
