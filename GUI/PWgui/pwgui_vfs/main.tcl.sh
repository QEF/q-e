cat > pwgui.vfs/main.tcl <<EOF

# load a starkit ...

package require starkit 
starkit::startup 

# load a Tk and Itcl ...

package require Tk $TK_VERSION      
package require $ITCL_EXACT Itcl $ITCL_VERSION

# manage the PWgui ...

puts " =================================================="
puts "  This is PWgui version: $PWGUI_VERSION"
puts " --------------------------------------------------"
puts " "


set pwgui app-pwgui
set guib  Guib-$GUIB_VERSION
set env(PWGUI) [file join \$starkit::topdir lib \$pwgui]
set env(GUIB)  [file join \$starkit::topdir lib \$pwgui lib \$guib]
source [file join \$starkit::topdir lib \$pwgui pwgui.tcl]
EOF
