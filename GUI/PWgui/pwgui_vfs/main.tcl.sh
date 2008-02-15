cat > $VFSDIR/main.tcl <<EOF

# load a starkit ...

package require starkit 
starkit::startup 

# load a Tk and Itcl ...

package require Tk $TK_VERSION      
package require $ITCL_EXACT Itcl $ITCL_VERSION

# manage the PWgui ...

puts " =================================================="
puts "  This is PWgui version: 3.1"
puts " --------------------------------------------------"
puts " "


set pwgui app-pwgui
set guib  Guib-0.3.4
set env(PWGUI) [file join \$starkit::topdir lib \$pwgui]
set env(GUIB)  [file join \$starkit::topdir lib \$guib]
source [file join \$starkit::topdir lib \$pwgui pwgui.tcl]
EOF
