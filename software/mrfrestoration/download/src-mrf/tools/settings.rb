# ********************************************************************
# qtsettings.rb
# A global file with the compiler and qt settings
#
# Author: Christian Wolf, christian.wolf@insa-lyon.fr
# ********************************************************************

$COMPILER="gcc"
$MAKE="make"

$MOC="moc"
$QTDIR="/usr/lib/qt4"
$QTLIBS="-L"+$QTDIR+"/lib -L/usr/X11R6/lib -lqt-mt -lX11 -lXext"
